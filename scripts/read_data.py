''' Reads wiggle files with ribosome profiling coverage data (from Shoelaces).
	Data can be in single file or split by length (and saved in a 'fwd' and 'rev' directories.)
	Maps the coverage to transcripts annotations (from read_genome_data.py).
'''

# fix 0/1-based

import numpy as np
from bx.binned_array import BinnedArray
import os
# import _pickle as pickle
# import gzip


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
	wig.close()
	return ribo_cov



def extract_intervals_from_wigs(wig_fwd, wig_rev, tx_coord, exp_exons=False):
	""" Dumps transcript coverage into python dictionary.
		tx_coord are transcript annotations from read_GTF/read_BED function.
		exp_exons=True returns full tx coverage, exons=False returns CDS coverage. """
 
	wf = read_WIG(wig_fwd)
	wr = read_WIG(wig_rev)
 
	orf_cov = {}
	# for every transcript (from BED)
	for tx in tx_coord:
		# strand = genomic_coord[tr][0].strand
		strand = tx_coord[tx]['strand']
		if strand == '+':
			handle = wf
		elif strand == '-':
			handle = wr
 
		# for every exon on transcript
		exons = np.ndarray(0)
		if not tx_coord[tx]['chromosome'] in handle.keys():
			continue
		else:
			chrom = tx_coord[tx]['chromosome']
 
		for i in range(len(tx_coord[tx]['exon_starts'])):
			start = tx_coord[tx]['exon_starts'][i]
			end = tx_coord[tx]['exon_ends'][i]
			exonint = handle[chrom].get_range(start, end+1) # +/-1 ??
			exons = np.append(exons, exonint)
 
		# if all values are nan, ignore transcript
		if all(np.isnan(exons)):
			continue
 
		# flip for '-' strand
		if strand == '-':
			exons = exons[::-1]
 
		### get ORF positions
		cds = exons[tx_coord[tx]['cds_coord'][0]:tx_coord[tx]['cds_coord'][1]]
 
		# if all values are nan, ignore transcript
		if all(np.isnan(cds)):
			continue
 
		### put all into a dictionary, return CDS coverage by default
		if exp_exons == True:
			orf_cov[tx] = exons
		else:
			orf_cov[tx] = cds
 
	return orf_cov


def extract_intervals_from_wigs_per_length(dirpath, tx_coord, exp_exons=True):
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
	for tx in tx_coord:
		# strand = genomic_coord[tr][0].strand
		strand = tx_coord[tx]['strand']
		if strand == '+':
			handles = fwd_handles
		elif strand == '-':
			handles = rev_handles
 
		exons_per_len = np.ndarray(0)
		if not tx_coord[tx]['chromosome'] in handle.keys():
			continue
		else:
			chrom = tx_coord[tx]['chromosome']
 
		# for different read lengths
		for i in range(len(handles)):
			exons = np.ndarray(0)
			for i in range(len(tx_coord[tx]['exon_starts'])):
				start = tx_coord[tx]['exon_starts'][i]
				end = tx_coord[tx]['exon_ends'][i]
				exonint = handle[chrom].get_range(start, end+1) # +/-1 ??
				exons = np.append(exons, exonint)
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
		cds = exons[:, tx_coord[tx]['cds_coord'][0]:tx_coord[tx]['cds_coord'][1]]
 
		# if all values are nan, ignore transcript
		if np.ndarray.all(np.isnan(cds)):
			continue
 
		### put all into a dictionary, return CDS coverage by default
		if exp_exons == True:
			orf_cov[tx] = exons_per_len
		else:
			orf_cov[tx] = cds
 
	return orf_cov



#######
# WIG - one file

wig_fwd = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-forward.wig'
wig_rev = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-reverse.wig'

wig_cov = extract_intervals_from_wigs(wig_fwd, wig_rev, transcripts)


