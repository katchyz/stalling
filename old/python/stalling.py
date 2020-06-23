import os
import cPickle as pickle
import gzip
import numpy as np
from scipy import sparse, stats
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals

# BEDfile = ...

def get_well_expressed_transcripts(data_dict, list_of_transcripts):
	""" Given the dictionary of all stages and transcripts, and a list of transcripts,
	(list of transcripts: orfs from shoelaces_full_distribution)
	returns the dictionary of highly expressed transcripts per each stage """
	codon_cov_high = {}
	for stage in data_dict:
		codon_cov_high[stage] = {}
		for tr in list_of_transcripts:
			nt_cov = np.asarray(data_dict[stage][tr][1].sum(axis=0))[0]
			# calculate codon coverage
			if (len(nt_cov) % 3 == 0) and (len(nt_cov) > 33):
				pep_array = np.zeros(len(nt_cov)/3)
				for i in range(0, len(nt_cov), 3):
					pep_array[i/3] = nt_cov[i] + nt_cov[i+1] + nt_cov[i+2]
				# choose CDSs where at least 50% of codons have at least one footprint
				if np.median(pep_array) > 0:
					codon_cov_high[stage][tr] = pep_array
	return codon_cov_high


def calculate_zscores(codon_cov_high):
	""" Given the dictionary of highly expressed transcripts, calculates
	the z-scores, avoiding first and last codons """
	zscores = {}
	for stage in codon_cov_high:
		zscores[stage] = {}
		for tr in codon_cov_high[stage]:
			zs = stats.zscore(codon_cov_high[stage][tr][1:-2])
			zs = np.insert(zs, 0, 0) # add a zero for 1st codon
			zs = np.append(zs, [0, 0]) # add a zero for last two codons
			zscores[stage][tr] = zs
	return zscores

def get_genes_with_peaks(zscores, threshold):
	""" Given the dictionary with z-scores and a threshold, return the dictionary
	of transcripts with list of peaks, for each stage separately """
	genes_with_peaks = {}
	for stage in zscores:
		genes_with_peaks[stage] = {}
		for tr in zscores[stage]:
			if max(zscores[stage][tr]) > threshold:
				peak_positions = []
				for i in range(len(zscores[stage][tr])):
					if zscores[stage][tr][i] > threshold:
						peak_positions += [i*3] # 1st nt in codon
				genes_with_peaks[stage][tr] = peak_positions
	return genes_with_peaks


def get_FASTA_sequence(filepath):
	""" Reads gzipped fasta file """
	f = gzip.open(filepath, 'rb')
	fasta = {}
	header_cds = ''
	seq_cds = ''
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if header_cds and seq_cds:
				fasta[header_cds] = seq_cds
			header_cds = line.split()[0][1:]
			seq_cds = ''
		else:
			seq_cds += line
		# for the last sequence:
		fasta[header_cds] = seq_cds
	return fasta


def get_FASTA_pep_sequence(filepath):
	""" Reads gzipped fasta file """
	f = gzip.open(filepath, 'rb')
	peptide = {}
	header_pep = ''
	seq_pep = ''
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if header_pep and seq_pep:
				peptide[header_pep] = seq_pep
			header_pep = line.split()[4][11:]
			seq_pep = ''
		else:
			seq_pep += line
		# for the last sequence:
		peptide[header_pep] = seq_pep
	return peptide


def get_fragment_lengths_around_stall_sites(data_dict, genes_with_peaks, no_lengths):
	""" Given dictionary with transcript data and set of genes with peaks,
	returns an array of fragments of different lengths around stall sites """
	meta_stall = {}
	for stage in genes_with_peaks:
		meta_stall[stage] = np.zeros((no_lengths,33)) ### for 7 different lengths
		for tr in genes_with_peaks[stage]:
			for peak in genes_with_peaks[stage][tr]:
				start = peak-15
				end = peak+18
				for l in range(no_lengths):
					if start < 0:
						meta_stall[stage][l][abs(start):] += data_dict[stage][tr][1].toarray()[l][0:end]
					elif end > len(data_dict[stage][tr][1].toarray()[l]):
						meta_stall[stage][l][0:33-(end-len(data_dict[stage][tr][1].toarray()[l]))] += data_dict[stage][tr][1].toarray()[l][start:len(data_dict[stage][tr][1].toarray()[l])]
					else:
						meta_stall[stage][l] += data_dict[stage][tr][1].toarray()[l][start:end]
	return meta_stall


def get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, upstream, downstream):
	""" Gets fasta sequences and list of stall sites, and number of nucleotides or amino acids
	upstream and downstream of the stall codon, returns lists of sequences for each stage.
	Upstream should be a negative integer, downstream a positive """
	seq_logo = {}
	for stage in genes_with_peaks:
		seq_logo[stage] = []
		for tr in genes_with_peaks[stage]:
			if tr in fasta:
				for peak in genes_with_peaks[stage][tr]:
					start = peak + upstream ### vary the length
					end = peak + 3 + downstream ### vary the length
					if start < 0:
						seq_logo[stage].append('N'*abs(start)+fasta[tr][0:end])
					elif end > len(fasta[tr]):
						seq_logo[stage].append(fasta[tr][start:len(fasta[tr])]+'N'*(end-len(fasta[tr])))
					else:
						seq_logo[stage].append(fasta[tr][start:end])
	return seq_logo


# write peaks to BED file
def write_peaks_to_BED_file(genes_with_peaks, filepath, geneToInterval, geneLength, orfs):
	""" Given genes with peaks and output file, writes a BED file for different stages """
	bed = open(filepath, 'w')
	bed.write('track name="Stall sites" visibility=2 itenRgb="On"\n')
	for stage in genes_with_peaks:
		for tr in genes_with_peaks[stage]:
			# stitch the geneToInterval intervals into mRNA with genomic positions
			mRNA = []
			for exon in range(len(geneToInterval[tr])):
				mRNA.extend(range(geneToInterval[tr][exon].start, geneToInterval[tr][exon].end))
			# offset by CDS start position
			if geneToInterval[tr][0].strand == '+':
				CDS_offset = orfs[tr][0][0]
			elif geneToInterval[tr][0].strand == '-':
				CDS_offset = geneLength[tr] - orfs[tr][0][1] + 3
			for i in range(len(genes_with_peaks[stage][tr])):
				# offset by genes_with_peaks_uniq[stage][tr][i]
				if geneToInterval[tr][0].strand == '+':
					peak_pos = CDS_offset + genes_with_peaks[stage][tr][i]
					# get genomic positions
					gen_st_pos = mRNA[peak_pos]
					gen_end_pos = mRNA[peak_pos]+3 # w/o checking for split over exons...
					# print in bed
					bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'200,100,0'+'\n')
				elif geneToInterval[tr][0].strand == '-':
					peak_pos = CDS_offset + genes_with_peaks[stage][tr][i] - 3
					gen_st_pos = mRNA[-peak_pos]-3
					gen_end_pos = mRNA[-peak_pos]
					bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'0,100,200'+'\n')







