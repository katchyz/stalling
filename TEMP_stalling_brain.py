import os
import cPickle as pickle
import gzip
import numpy as np
from scipy import sparse, stats
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
from stalling import *

data_dict = seq_orf
#data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/Brain.p.gz', 'rb'))

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

codon_cov_high = get_well_expressed_transcripts(data_dict, orfs)

zscores = calculate_zscores(codon_cov_high)

genes_with_peaks = get_genes_with_peaks(zscores, 8.0)

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz')

meta_stall = get_fragment_lengths_around_stall_sites(data_dict, genes_with_peaks, 7)

seq_logo = get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, -15, 18)



################

#f = open('/Volumes/USELESS/META/LOGO_brain.txt', 'w')

for stage in seq_logo:
	f = open('/Volumes/USELESS/META/LOGO_Giraldez/'+stage+'.txt', 'w')
	for seq in seq_logo[stage]:
		if len(seq) == 36:
			f.write(seq+'\n')
	f.close()



#f = open('/Users/kasia/Downloads/mouse_gtp.txt', 'r')
#tg transcript: gene
gwp = {}
gwp['eb_hiseq'] = {}
genes = []
for tr in genes_with_peaks['eb_hiseq']:
	if tr in tg:
		gene = tg[tr]
		if gene not in genes:
			gwp['eb_hiseq'][tr] = genes_with_peaks['eb_hiseq'][tr]
			genes.append(gene)





for stage in seq_logo:
	f = open('/Volumes/USELESS/META/LOGO_yeast/'+stage+'.txt', 'w')
	for seq in seq_logo[stage]:
		if len(seq) == 36:
			f.write(seq+'\n')
	f.close()








