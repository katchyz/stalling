### FRUIT FLY ###

from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord
import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Drosophila_melanogaster.BDGP6.79.bed'
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



stage_dir_names = ['0-2h_embryo_A', '0-2h_embryo_B', '0-2h_embryo_cushion_A', '0-2h_embryo_cushion_B', \
'0-2h_embryo_gradient_A', '0-2h_embryo_gradient_B', 'S2Cell_150_A', 'S2cell_150_B', 'S2cell_250_A', 'S2cell_250_B']

stage_dir_28_37 = ['0-2h_embryo_A', '0-2h_embryo_B']
stage_dir_28_36 = ['0-2h_embryo_cushion_A', '0-2h_embryo_cushion_B', '0-2h_embryo_gradient_A', '0-2h_embryo_gradient_B']
stage_dir_27_37 = ['S2Cell_150_A', 'S2cell_150_B']
stage_dir_27_35 = ['S2cell_250_A', 'S2cell_250_B']

for stagedir in stage_dir_27_35:
	for fr_len in range(27, 36):
		main_dir = '/Volumes/USELESS/OUT/fruitfly_out/by_length/'+stagedir
		tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval, stagedir, fr_len)
		pickle.dump(tr_len, gzip.open('/Volumes/USELESS/META/fruitfly_meta/'+stagedir+'-'+str(fr_len)+'.p.gz', 'wb'))
		print 'saved ', stagedir, fr_len, '\n'
		tr_len = None


