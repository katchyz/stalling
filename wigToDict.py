from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord
import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip

bedpath = '/Users/kasia/Documents/PhD/data/Danio_rerio.Zv9.79_chr.bed'
bedpath = '/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)


def extract_intervals_from_wigs(dirpath, list_of_transcripts, geneToInterval):
	""" Dumps transcript coverage into python dictionary """
	fwd_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'fwd')))):
		fwd_handles.extend([ParseWig(os.path.join(dirpath, 'fwd', file))])
	rev_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'rev')))):
		rev_handles.extend([ParseWig(os.path.join(dirpath, 'rev', file))])
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
				exonint = np.asarray(handles[i].fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
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


#stage_dir_names = ['0_2-4Cell', '1_256Cell', '2_1KCell', '3_Dome', '4_Shield', '5_Bud', '6_28hpf', '7_5dpf']
stage_dir_names = ['0_2-4Cell', '1_256Cell', '2_1KCell', '5_Bud']

for stagedir in stage_dir_names:
	main_dir = '/Users/kasia/Documents/PhD/outputs/Zv9/our_remapped/aln3/by_length/'+stagedir
	tr_len = extract_intervals_from_wigs(main_dir, orfs, geneToInterval)
	pickle.dump(tr_len, gzip.open('/Users/kasia/Documents/PhD/meta_data/Zv9/'+stagedir+'.p.gz', 'wb'))
	print 'saved ', stagedir
	tr_len = None


############################
############################
############################
############################


from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord
import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip

bedpath = '/Users/kasia/Documents/PhD/data/Danio_rerio.Zv9.79_chr.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

def extract_intervals_from_wigs_fwd(dirpath, list_of_transcripts, geneToInterval):
	""" Dumps transcript coverage into python dictionary """
	fwd_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'fwd')))):
		fwd_handles.extend([ParseWig(os.path.join(dirpath, 'fwd', file))])
	transcripts_lengths = {}
	for tr in list_of_transcripts:
		strand = geneToInterval[tr][0].strand
		if strand == '+':
			handles = fwd_handles
			sequencearray = np.zeros((len(fwd_handles), geneLength[tr])) # IF USING FWD_HANDLES!
			for i in range(len(handles)): # for different read lengths
				cds = np.array([])
				for interval in geneToInterval[tr]:
					exonint = np.asarray(handles[i].fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
					cds = np.append(cds, exonint)
				cds[np.isnan(cds)] = 0 # get all exons for one handle
				sequencearray[i] = cds ###lengths???
			### get ORF positions
			orfarray = sequencearray[:, orfs[tr][0][0]:orfs[tr][0][1]]
		# for + strand, sparse
			sp_orf = sparse.csr_matrix(orfarray)
			sp_tran = sparse.csr_matrix(sequencearray)
			transcripts_lengths[tr] = sp_tran, sp_orf
	return transcripts_lengths

def extract_intervals_from_wigs_rev(dirpath, list_of_transcripts, geneToInterval):
	""" Dumps transcript coverage into python dictionary """
	rev_handles = []
	for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(dirpath, 'rev')))):
		rev_handles.extend([ParseWig(os.path.join(dirpath, 'rev', file))])
	transcripts_lengths = {}
	for tr in list_of_transcripts:
		strand = geneToInterval[tr][0].strand
		if strand == '-':
			handles = rev_handles
		### initiate sequence array
			sequencearray = np.zeros((len(rev_handles), geneLength[tr])) # IF USING FWD_HANDLES!
			for i in range(len(handles)): # for different read lengths
				cds = np.array([])
				for interval in geneToInterval[tr]:
					exonint = np.asarray(handles[i].fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
					cds = np.append(cds, exonint)
				cds[np.isnan(cds)] = 0 # get all exons for one handle
				sequencearray[i] = cds ###lengths???
			### get ORF positions
			orfarray = sequencearray[:, orfs[tr][0][0]:orfs[tr][0][1]]
			orfarray = np.fliplr(orfarray)
			sequencearray = np.fliplr(sequencearray)
			sp_orf = sparse.csr_matrix(orfarray)
			sp_tran = sparse.csr_matrix(sequencearray)
		### put all into a dictionary
			transcripts_lengths[tr] = sp_tran, sp_orf
	return transcripts_lengths

##################
##################

stage_dir_names = ['3_Dome', '4_Shield', '6_28hpf', '7_5dpf']

for stagedir in stage_dir_names:
	main_dir = '/Users/kasia/Documents/PhD/outputs/Zv9/our_remapped/aln3/by_length/'+stagedir
	tr_len = extract_intervals_from_wigs_fwd(main_dir, orfs, geneToInterval)
	pickle.dump(tr_len, gzip.open('/Users/kasia/Documents/PhD/meta_data/Zv9/'+stagedir+'_fwd.p.gz', 'wb'))
	print 'saved ', stagedir
	tr_len = None


stage_dir_names = ['3_Dome', '4_Shield', '6_28hpf', '7_5dpf']

for stagedir in stage_dir_names:
	main_dir = '/Users/kasia/Documents/PhD/outputs/Zv9/our_remapped/aln3/by_length/'+stagedir
	tr_len = extract_intervals_from_wigs_rev(main_dir, orfs, geneToInterval)
	pickle.dump(tr_len, gzip.open('/Users/kasia/Documents/PhD/meta_data/Zv9/'+stagedir+'_rev.p.gz', 'wb'))
	print 'saved ', stagedir
	tr_len = None


