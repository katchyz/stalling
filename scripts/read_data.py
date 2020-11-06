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
 
	geneToInterval = {}
	cds = {}
 
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
 
		if geneStr in geneToInterval:
			sys.stderr.write("Warning: duplicate gene names: %s. Overwriting previous entry.\n" % geneStr)
		geneToInterval[geneStr] = []
 
		offset = 0
		for exonStart, exonLen in zip(exonStarts, exonLengths):
			exonStart = int(exonStart)
			exonLen = int(exonLen)
			s = tx_start + exonStart
			e = s + exonLen
 
			# genomic intervals
			interval = GenomicInterval(chrom, s, e, strand)
			geneToInterval[geneStr].append(interval)
 
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
			cds[geneStr] = []
		cds[geneStr].append((new_cds_start, new_cds_end))
 
	return geneToInterval, cds



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


def fetch_scores(ribo_cov, chr, st, end):
	if chr in ribo_cov.keys():
		return [ ribo_cov[chr][i] for i in range(st,end)]
	else:
		return [np.nan]*(end-st)


def extract_intervals_from_wigs_per_length(dirpath, list_of_transcripts, geneToInterval):
	""" Dumps transcript coverage into python dictionary.
		dirpath is a path to directory with two directories: 'fwd' and 'rev' containing wig files split by length.
		list_of_transcripts (cds) and geneToInterval are from read_BED function. """
	fwd_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'fwd')))):
		fwd_handles.extend([read_WIG(os.path.join(dirpath, 'fwd', file))])
	rev_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'rev')))):
		rev_handles.extend([read_WIG(os.path.join(dirpath, 'rev', file))])
	transcripts_lengths = {}
	for tr in list_of_transcripts:
		strand = geneToInterval[tr][0].strand
		if strand == '+':
			handles = fwd_handles
		elif strand == '-':
			handles = rev_handles
		### initiate sequence array
		sequencearray = np.zeros((len(fwd_handles), geneLength[tr])) # IF USING FWD_HANDLES!
		for i in range(len(handles)): # for different read lengths
			cds = np.array([])
			for interval in geneToInterval[tr]:
				exonint = np.asarray(fetch_scores(handles[i], interval.chrom, interval.start+1, interval.end+1))
				cds = np.append(cds, exonint)
			cds[np.isnan(cds)] = 0 # get all exons for one handle
			sequencearray[i] = cds ###lengths???
		### get ORF positions
		orfarray = sequencearray[:, orfs[tr][0][0]:orfs[tr][0][1]]
		# for + strand, sparse
		sp_orf = sparse.csr_matrix(orfarray)
		sp_tran = sparse.csr_matrix(sequencearray)
		if strand == '-':
			orfarray = np.fliplr(orfarray)
			sequencearray = np.fliplr(sequencearray)
			sp_orf = sparse.csr_matrix(orfarray)
			sp_tran = sparse.csr_matrix(sequencearray)
		### put all into a dictionary
		transcripts_lengths[tr] = sp_tran, sp_orf
	return transcripts_lengths


bed = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
wig_fwd = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-forward.wig'
wig_rev = '/Volumes/USELESS/STALLING/wigs/one_file/mouse_ingolia2011_NONE-reverse.wig'

(geneToInterval, orfs) = BED_to_intervals(bed)
main_dir = '/Volumes/USELESS/outputs_MAC/Zv9/our_remapped/aln3/by_length/5_Bud'
tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval)


# speed up extract_intervals from wigs
# make options for single and split by length wig files
