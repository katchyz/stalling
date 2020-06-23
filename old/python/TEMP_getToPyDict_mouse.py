### MOUSE ###

from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord
import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)



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
		seqarray = np.array([])
		for interval in geneToInterval[tr]:
			exonint = np.asarray(handle.fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
			seqarray = np.append(seqarray, exonint)
		seqarray[np.isnan(seqarray)] = 0 # get all exons for one handle
		### get ORF positions
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



stage_dir_names = ['EB_GA_1', 'EB_GA_2', 'EB_GA_3', 'EB_HiSeq', 'ESC_36h_GA_1', 'ESC_36h_GA_2', 'ESC_36h_GA_3', \
'ESC_36h_HiSeq', 'ESC_ff_GA_1', 'ESC_ff_GA_2', 'ESC_ff_GA_3', 'ESC_ff_HiSeq', 'ESC_ff_chx_GA_1', 'ESC_ff_chx_GA_2', \
'ESC_ff_emet_GA', 'ESC_ff_none_GA_1', 'ESC_ff_none_GA_2', 'ESC_ff_none_GA_3', 'ESC_ff_none_GA_4']

stage_dir_names = ['ESC_ff_GA_3', 'ESC_ff_HiSeq', 'ESC_ff_chx_GA_1', 'ESC_ff_chx_GA_2', \
'ESC_ff_emet_GA', 'ESC_ff_none_GA_1', 'ESC_ff_none_GA_2', 'ESC_ff_none_GA_3', 'ESC_ff_none_GA_4']

stage_dir_names = ['ESC_ff_none_GA_1', 'ESC_ff_none_GA_2', 'ESC_ff_none_GA_3', 'ESC_ff_none_GA_4']

stage_dir_names = ['ESC_ff_GA_3', 'ESC_ff_HiSeq', 'ESC_ff_chx_GA_1', 'ESC_ff_chx_GA_2', 'ESC_ff_emet_GA']

stage_dir_names = ['mouse_none']

stage_dir_names = ['all_len_no_off']

for stagedir in stage_dir_names:
	for fr_len in range(25, 36):
		main_dir = '/Volumes/USELESS/OUT/mouse_out/by_length/'+stagedir
		tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval, stagedir, fr_len)
		pickle.dump(tr_len, gzip.open('/Volumes/USELESS/META/mouse_meta/'+stagedir+'-'+str(fr_len)+'.p.gz', 'wb'))
		print 'saved ', stagedir, fr_len, '\n'
		tr_len = None


stage_dir_names = ['Brain1_pm', 'Brain2_mn', 'Brain3_am', 'Brain4_nn']

for stagedir in stage_dir_names:
	for fr_len in range(25, 32):
		main_dir = '/Volumes/USELESS/OUT/zebrafish_brain_out/by_length/'+stagedir
		tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval, stagedir, fr_len)
		pickle.dump(tr_len, gzip.open('/Volumes/USELESS/META/zebrafish_brain_meta/'+stagedir+'-'+str(fr_len)+'.p.gz', 'wb'))
		print 'saved ', stagedir, fr_len, '\n'
		tr_len = None


stage_dir_names = ['mouse_all_len_no_off']

for stagedir in stage_dir_names:
	for fr_len in [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42]:
		main_dir = '/Volumes/USELESS/OUT/mouse_out/all_len_no_off'
		tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval, stagedir, fr_len)
		pickle.dump(tr_len, gzip.open('/Volumes/USELESS/META/mouse_meta/mouse'+stagedir+'-'+str(fr_len)+'.p.gz', 'wb'))
		print 'saved ', stagedir, fr_len, '\n'
		tr_len = None



#############

import numpy as np
import cPickle as pickle
import gzip

d = {}

stage_dir_names = ['mouse_all_len_no_off']

for stage in stage_dir_names:
	d[stage] = {}
	for l in [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42]:
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/mouse_meta/mouse'+stage+'-'+str(l)+'.p.gz', 'rb'))

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
	for tr in d[stage]['20']:
		seq_orf[stage][tr] = d[stage]['20'][tr][0], d[stage]['20'][tr][1]

for stage in seq_orf:
	for l in [21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42]:
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])


pickle.dump(seq_orf, gzip.open('/Volumes/USELESS/META/mouse_all_len_no_off.p.gz', 'wb'))
