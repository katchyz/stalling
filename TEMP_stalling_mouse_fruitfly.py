# MOUSE - none - stalling

import os
import cPickle as pickle
import gzip
import numpy as np
from scipy import sparse, stats
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
from stalling import *

data_dict = {}
data_dict['ESC_none'] = ESC_none
#data_dict['ESC_none'] = pickle.load(gzip.open('/Volumes/USELESS/META/ESC_none_added.p.gz', 'rb'))

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Drosophila_melanogaster.BDGP6.79.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

codon_cov_high = get_well_expressed_transcripts(data_dict, orfs)

zscores = calculate_zscores(codon_cov_high)

genes_with_peaks = get_genes_with_peaks(zscores, 8.0)

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')
fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz')

meta_stall = get_fragment_lengths_around_stall_sites(data_dict, genes_with_peaks)

seq_logo = get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, -15, 18)

seq_logo = get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, -21, 21)

################

f = open('/Volumes/USELESS/META/LOGO_ESC_none_-21+24.txt', 'w')
for seq in seq_logo['ESC_none']:
	if len(seq) == 45:
		f.write(seq+'\n')
f.close()

f = open('/Volumes/USELESS/META/LOGO_0-2h_embryo_A-21+24.txt', 'w')
for seq in seq_logo['0-2h_embryo_A']:
	if len(seq) == 45:
		f.write(seq+'\n')

f = open('/Volumes/USELESS/META/LOGO_S2cell_250_B-21+24.txt', 'w')
for seq in seq_logo['S2cell_250_B']:
	if len(seq) == 45:
		f.write(seq+'\n')

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














