# genes with peaks - zebrafish, mouse, fruitfly

import os
import cPickle as pickle
import gzip
import numpy as np
from scipy import sparse, stats
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
from stalling import *


data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/mouse_none.p.gz', 'rb'))
# data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/fruitfly_meta/fruitfly.p.gz', 'rb'))
# data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/Giraldez.p.gz', 'rb'))
# data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/OUR.p.gz', 'rb'))
# data_dict = pickle.load(gzip.open('/Volumes/USELESS/META/yeast.p.gz', 'rb'))

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
# bedpath = '/Volumes/USELESS/DATA/genomes/BED/Drosophila_melanogaster.BDGP6.79.bed'
# bedpath = '/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed'
#
# bedpath = '/Volumes/USELESS/DATA/genomes/BED/Saccharomyces_cerevisiae.R64-1-1.79.bed'

(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

codon_cov_high = get_well_expressed_transcripts(data_dict, orfs)

zscores = calculate_zscores(codon_cov_high)

genes_with_peaks = get_genes_with_peaks(zscores, 5.0)

# fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')
# meta_stall = get_fragment_lengths_around_stall_sites(data_dict, genes_with_peaks)
# seq_logo = get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, -15, 18)

pickle.dump(genes_with_peaks, gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_mouse.p.gz', 'wb'))
# pickle.dump(genes_with_peaks, gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_fruitfly.p.gz', 'wb'))
# pickle.dump(genes_with_peaks, gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_Giraldez.p.gz', 'wb'))
# pickle.dump(genes_with_peaks, gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_OUR.p.gz', 'wb'))
# pickle.dump(genes_with_peaks, gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_yeast.p.gz', 'wb'))

################

f = open('/Volumes/USELESS/META/genes_with_peaks/yeast.txt', 'w')
for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		f.write(tr+'\n')
f.close















