### Call peaks (putative stall sites) from ribosome coverage data

## read BED (GTF?) file with genomic coordinates
## read wig file with ribo coverage

from HTSeq import GenomicInterval
import os
import numpy as np
from scipy import sparse
from bx.binned_array import BinnedArray

# import _pickle as pickle
# import gzip


def read_BED(bed_file):
	''' Reads genomic and mRNA coordinates from BED file.
		Returns dictionary of genomic intervals (exons) for each transcript
		and dictionary of CDS start and stop on mRNA coordinates. '''
 
	bed = open(bed_file, "r")
 
	genomic_coord = {}
	cds_coord = {}
 
	for line in bed:
		gene = line.split("\t")
		chrom = gene[0]
		tx_start = int(gene[1])
		geneStr = gene[3] #"%s:::%s" % (gene[3], gene[1])
		strand = gene[5]
		cds_start = int(gene[6]) - tx_start
		cds_end = int(gene[7]) - tx_start
 
		# Non-coding
		if cds_start == cds_end:
			continue
		
		exonLengths = gene[10].split(",")
		exonStarts = gene[11].split(",")
 
		# BED files occasionally end with a ","
		if gene[10][-1] == ",":
			exonLengths.pop()
		if gene[11][-1] == ",":
			exonStarts.pop()
 
		if geneStr in genomic_coord:
			sys.stderr.write("Warning: duplicate gene names: %s. Overwriting previous entry.\n" % geneStr)
		genomic_coord[geneStr] = []
 
		offset = 0
		for exonStart, exonLen in zip(exonStarts, exonLengths):
			exonStart = int(exonStart)
			exonLen = int(exonLen)
			s = tx_start + exonStart
			e = s + exonLen
 
			# genomic intervals
			interval = GenomicInterval(chrom, s, e, strand)
			genomic_coord[geneStr].append(interval)
 
			# mRNA coordinates
			if cds_start >= exonStart and cds_start <= exonStart+exonLen:
				new_cds_start = offset + (cds_start - exonStart)
 
			if cds_end >= exonStart and cds_end <= exonStart+exonLen:
				new_cds_end = offset + (cds_end - exonStart)
 
			offset += exonLen
 
		if new_cds_start > new_cds_end:
			sys.stderr.write("CDS length is negative for transcript %s\n" % gene[3])
			continue
 
		if not geneStr in cds:
			cds_coord[geneStr] = []
		cds_coord[geneStr].append((new_cds_start, new_cds_end))
 
	return genomic_coord, cds_coord



def read_WIG(wig_file):
	wig = open(wig_file, "r")
	ribo_cov = {}
	for line in wig:
		if line.startswith('variableStep'):  # specific for Shoelaces wig files
			stepType = 'variable'
			fields = line.split()[1:]
			declarations = dict([(p[0], p[1].strip('"')) for p in [x.split("=") for x in fields]])
			chrom = declarations['chrom']
			span = 1
			ribo_cov[chrom] = BinnedArray()
		else:
			tmp = line.strip().split()
			pos = int(tmp[0])
			val = float(tmp[1])
			ribo_cov[chrom][pos] = val
	return ribo_cov


# def fetch_scores(ribo_cov, chr, st, end):
# 	if chr in ribo_cov.keys():
# 		return [ ribo_cov[chr][i] for i in range(st,end)]
# 	else:
# 		return [np.nan]*(end-st)


def extract_intervals_from_wigs(wig_fwd, wig_rev, cds_coord, genomic_coord, exp_exons=False):
	""" Dumps transcript coverage into python dictionary.
		dirpath is a path to directory with two directories: 'fwd' and 'rev' containing wig files split by length.
		CDS_coord and genomic_coord are from read_BED function.
		exp_exons=True returns full tx coverage, exons=False returns CDS coverage. """
 
	wf = read_WIG(wig_fwd)
	wr = read_WIG(wig_rev)
 
	orf_cov = {}
	# for every transcript (from BED)
	for tr in cds_coord:
		strand = genomic_coord[tr][0].strand
		if strand == '+':
			handle = wf
		elif strand == '-':
			handle = wr
 
		# for every exon on transcript
		exons = np.ndarray(0)
		for interval in genomic_coord[tr]:
			if interval.chrom in handle.keys():
				exonint = handle[interval.chrom].get_range(interval.start+1, interval.end+1)
				exons = np.append(exons, exonint)
				# cds[np.isnan(cds)] = 0 # get all exons for one handle
 
		# if all values are nan, ignore transcript
		if all(np.isnan(exons)):
			continue
 
		# flip for '-' strand
		if strand == '-':
			exons = exons[::-1]
 
		### get ORF positions
		cds = exons[cds_coord[tr][0][0]:cds_coord[tr][0][1]]
 
		# if all values are nan, ignore transcript
		if all(np.isnan(cds)):
			continue
 
		### put all into a dictionary, return CDS coverage by default
		if exp_exons == True:
			orf_cov[tr] = exons
		else:
			orf_cov[tr] = cds
 
	return orf_cov


def extract_intervals_from_wigs_per_length(dirpath, cds_coord, genomic_coord, exp_exons=True):
	""" Dumps transcript coverage into python dictionary.
		dirpath is a path to directory with two directories: 'fwd' and 'rev' containing wig files split by length.
		list_of_transcripts (cds) and geneToInterval are from read_BED function.
		exp_exons=True returns full tx coverage, exons=False returns CDS coverage. """
 
	fwd_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'fwd')))):
		fwd_handles.extend([read_WIG(os.path.join(dirpath, 'fwd', file))])
	rev_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'rev')))):
		rev_handles.extend([read_WIG(os.path.join(dirpath, 'rev', file))])
 
	orf_cov = {}
	for tr in cds_coord:
		strand = genomic_coord[tr][0].strand
		if strand == '+':
			handles = fwd_handles
		elif strand == '-':
			handles = rev_handles
 
		exons_per_len = np.ndarray(0)
		# for different read lengths
		for i in range(len(handles)):
			exons = np.ndarray(0)
			for interval in genomic_coord[tr]:
				if interval.chrom in handles[i].keys():
					exonint = handles[i][interval.chrom].get_range(interval.start+1, interval.end+1)
					exons = np.append(exons, exonint)
					# exons[np.isnan(exons)] = 0 # get all exons for one handle
			if i == 0:
				exons_per_len = np.append(exons_per_len, exons)
			else:
				np.vstack((exons_per_len, exons))
 
		# if all values are nan, ignore transcript
		if np.ndarray.all(np.isnan(exons_per_len)):
			continue
 
		# flip for '-' strand
		if strand == '-':
			exons_per_len = np.flip(exons_per_len)
 
		### get ORF positions
		cds = exons[:, cds_coord[tr][0][0]:cds_coord[tr][0][1]]
 
		# if all values are nan, ignore transcript
		if np.ndarray.all(np.isnan(cds)):
			continue
 
		### put all into a dictionary, return CDS coverage by default
		if exp_exons == True:
			orf_cov[tr] = exons_per_len
		else:
			orf_cov[tr] = cds
 
	return orf_cov





bed = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
wig_fwd = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-forward.wig'
wig_rev = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-reverse.wig'

(genomic_coord, cds_coord) = read_BED(bed)
wig_cov = extract_intervals_from_wigs(wig_fwd, wig_rev, cds_coord, genomic_coord)


