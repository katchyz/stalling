# zebrafish_Giraldez

from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord
import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

stage_dir_names = ['1KCell', '2-4Cell', '256Cell', '28hpf', '5dpf', 'Bud', 'Dome', 'Shield', 'Shield_MZOep', 'Shield_Squint']


### BIG files: 'RPF_02h_run2', 'RPF_24h_run2'

def extract_intervals_from_wigs(dirpath, list_of_transcripts, geneToInterval, filename, fragment_length):
	""" Dumps transcript coverage into python dictionary """
	fl = str(fragment_length)
	fwd_handle = ParseWig(os.path.join(dirpath, 'fwd', filename+'-'+fl+'-forward.wig'))
	rev_handle = ParseWig(os.path.join(dirpath, 'rev', filename+'-'+fl+'-reverse.wig'))
	transcripts_lengths = {}
	for tr in list_of_transcripts:
		strand = geneToInterval[tr][0].strand
		if strand == '+':
			handle = fwd_handle
		elif strand == '-':
			handle = rev_handle
		### initiate sequence array
		#sequencearray = np.zeros((1, geneLength[tr]))
		seqarray = np.array([])
		for interval in geneToInterval[tr]:
			exonint = np.asarray(handle.fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
			seqarray = np.append(seqarray, exonint)
		seqarray[np.isnan(seqarray)] = 0 # get all exons for one handle
		#sequencearray[0] = cds ###lengths???
		### get ORF positions
		#orfarray = sequencearray[:, orfs[tr][0][0]:orfs[tr][0][1]]
		orfarray = seqarray[orfs[tr][0][0]:orfs[tr][0][1]]
		# for + strand, sparse
		sp_orf = sparse.csr_matrix(orfarray)
		sp_seq = sparse.csr_matrix(seqarray)
		if strand == '-':
			orfarray = orfarray[::-1]
			seqarray = seqarray[::-1]
			sp_orf = sparse.csr_matrix(orfarray)
			sp_seq = sparse.csr_matrix(seqarray)
		### put all into a dictionary
		transcripts_lengths[tr] = sp_seq, sp_orf
	return transcripts_lengths


for stagedir in stage_dir_names:
	for fr_len in range(25, 32):
		main_dir = '/Volumes/USELESS/OUT/zebrafish_our_out/by_length/'+stagedir
		tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval, stagedir, fr_len)
		pickle.dump(tr_len, gzip.open('/Volumes/USELESS/META/zebrafish_our_meta/'+stagedir+'-'+str(fr_len)+'.p.gz', 'wb'))
		print 'saved ', stagedir, fr_len, '\n'
		tr_len = None



##############
## merge GZ files
##############


import numpy as np
import cPickle as pickle
import gzip

d = {}

for stage in stage_dir_names:
	d[stage] = {}
	for l in range(25, 32):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/zebrafish_our_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

def csr_vappend(a,b):
	""" Takes in 2 csr_matrices and appends the second one to the bottom of the first one. 
	Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
	the first matrix instead of copying it. The data, indices, and indptr still get copied."""
	a.data = np.hstack((a.data,b.data))
	a.indices = np.hstack((a.indices,b.indices))
	a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
	a._shape = (a.shape[0]+b.shape[0],b.shape[1])
	return a


#

seq_orf = {}
for stage in d:
	seq_orf[stage] = {}
	for tr in d[stage]['25']:
		seq_orf[stage][tr] = d[stage]['25'][tr][0], d[stage]['25'][tr][1]

for stage in seq_orf:
	for l in range(26, 32):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])


pickle.dump(seq_orf, gzip.open('/Volumes/USELESS/META/OUR.p.gz', 'wb'))