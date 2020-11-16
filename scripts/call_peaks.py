''' Finds the longest transcript per gene (optionally),
	selects transcripts with good coverage (median codon coverage > 0),
	calculates z-score coverage,
	calls peaks for each transcript.
'''

import numpy as np
from scipy import stats



def get_longest_transcript(transcripts):
	""" Takes output from read_BED or read_GTF and returns only the longest tx (cds) per gene """
 
	longest_tx = {}
 
	for tx in transcripts:
		gene = transcripts[tx]['gene_id']
		cds_len = transcripts[tx]['cds_coord'][1] - transcripts[tx]['cds_coord'][0]
 
		if not gene in longest_tx.keys():
			longest_tx[gene] = transcripts[tx]
			longest_tx[gene]['tx_id'] = tx
		else:
			if cds_len > longest_tx[gene]['cds_coord'][1] - longest_tx[gene]['cds_coord'][0]:
				longest_tx[gene] = transcripts[tx]
				longest_tx[gene]['tx_id'] = tx
 
	return longest_tx


def get_well_expressed_transcripts(ribo_cov):
	""" Takes a dictionary of ribosome coverage and returns a dictionary with well expressed transcripts
		(median codon coverage less than one). """
 
	well_expr_tx = {}
 
	for tx in ribo_cov:
		# nt_cov = np.asarray(data_dict[stage][tr][1].sum(axis=0))[0] per length
		nt_cov = np.nan_to_num(ribo_cov[tx])
		# calculate codon coverage
		if len(nt_cov) > 33: # len(nt_cov) % 3 == 0
			pep_array = np.zeros(int(np.floor(len(nt_cov)/3)))
			for i in range(0, int(np.floor(len(nt_cov)/3)*3), 3):
				pep_array[int(i/3)] = nt_cov[i] + nt_cov[i+1] + nt_cov[i+2]
 
			# choose CDSs where at least 50% of codons have at least one footprint
			if np.median(pep_array) > 0:
				well_expr_tx[tx] = pep_array
 
		else:
			print(f'Length of gene {tx} is less than 33 nt')
 
	return well_expr_tx


def calculate_zscores(well_expr_tx, start=1, end=2):
	""" Given the dictionary of highly expressed transcripts, calculates
		the z-scores, avoiding first and last codons.
		start and end are numbers of codons at the beginning/end of tx to exclude from calculation. """
 
	zscores = {}
 
	for tx in well_expr_tx:
		zs = stats.zscore(well_expr_tx[tx][start:-end])
		zs = np.insert(zs, 0, [0]*start) # add a zero for 1st codon
		zs = np.append(zs, [0]*end) # add a zero for last two codons
		zscores[tx] = zs
 
	return zscores


def get_genes_with_peaks(zscores, threshold=5):
	""" Given the dictionary with z-scores and a z-score threshold,
		return the dictionary of transcripts with list of peaks. """
 
	genes_with_peaks = {}
 
	for tx in zscores:
		if max(zscores[tx]) > threshold:
			peak_positions = np.where(zscores[tx] > threshold)[0] *3 # 1st nt in codon
			genes_with_peaks[tx] = peak_positions
 
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





longest_tx = get_longest_transcript(transcripts)
ribo_cov = extract_intervals_from_wigs(wig_fwd, wig_rev, longest_tx)
wet = get_well_expressed_transcripts(ribo_cov)
zscores = calculate_zscores(wet, start=5, end=2)
gwp = get_genes_with_peaks(zscores, threshold=5)





