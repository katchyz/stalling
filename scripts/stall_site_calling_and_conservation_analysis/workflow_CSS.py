#!/usr/bin/python3

''' Workflow for finding CSSs '''

import os
import numpy as np
from bx.binned_array import BinnedArray
from scipy import stats
import gzip
import _pickle as pickle
import subprocess
from itertools import combinations
import pandas as pd

from read_genome_annotations import *
from read_data import *
from read_sequence import *
from call_peaks import *
from conserved_stall_sites import *
from save_data import *

### READ DATA

## GENOMIC COORDINATES (read_genome_annotations.py)
path_gtf = '../DATA/genomes/GTF'

gtf = {}
gtf['yeast'] = read_GTF(os.path.join(path_gtf, 'Saccharomyces_cerevisiae.R64-1-1.79.gtf'))
gtf['fruitfly'] = read_GTF(os.path.join(path_gtf, 'Drosophila_melanogaster.BDGP6.79.gtf'))
gtf['zebrafish'] = read_GTF(os.path.join(path_gtf, 'Danio_rerio.GRCz10.84.chr.gtf'))
gtf['mouse'] = read_GTF(os.path.join(path_gtf, 'Mus_musculus.GRCm38.79.chr.gtf'))
gtf['human'] = read_GTF(os.path.join(path_gtf, 'Homo_sapiens.GRCh38.79.chr.gtf'))

# get longest transcripts per gene
genes = {}
for key in gtf.keys():
	genes[key] = get_longest_transcript(gtf[key])

# change dictionary key from gene to tx
transcripts = {}
for key in genes.keys():
	transcripts[key] = set_keys_to_tx_id(genes[key])


## SEQUENCE DATA (read_sequence.py)
path = '../DATA'

fasta = {}
fasta['yeast'] = get_FASTA_sequence(os.path.join(path, 'fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz'))
fasta['fruitfly'] = get_FASTA_sequence(os.path.join(path, 'fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz'))
fasta['zebrafish'] = get_FASTA_sequence(os.path.join(path, 'fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz'))
# names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
for key in fasta['zebrafish'].keys():
	fasta['zebrafish'][key[0:18]] = fasta['zebrafish'][key]
	del fasta['zebrafish'][key]

fasta['mouse'] = get_FASTA_sequence(os.path.join(path, 'fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz'))
fasta['human'] = get_FASTA_sequence(os.path.join(path, 'fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz'))


## RIBOSOME PROFILING DATA (read_data.py)
path_wig = '../DATA/wigs/one_file'
wig_files = [f for f in os.listdir(path_wig) if os.path.isfile(os.path.join(path_wig, f))]
wig_unique = list(set([i.split('-')[0] for i in wig_files])) # unique prefixes of wig files
ribo_cov = {}
for prefix in wig_unique:
	org = prefix.split('_')[0]
	wig_cov = extract_intervals_from_wigs(os.path.join(cons_path, prefix, '-forward.wig'), \
		os.path.join(cons_path, prefix, '-reverse.wig'), transcripts[org])
	ribo_cov[prefix] = wig_cov

# get well-expressed transcripts
wet = {}
for lib in ribo_cov.keys():
	wet[lib] = get_well_expressed_transcripts(ribo_cov[lib])

# calculate z-scores
zscores = {}
for lib in wet.keys():
	zscores[lib] = calculate_zscores(wet[lib], start=5, end=2)

# get genes with peaks
gwp = {}
for lib in zscores.keys():
	gwp[lib] = get_genes_with_peaks(zscores, threshold=5)

# get peaks for each organism in separate library (save to CSV)
gwpeaks = {}
for lib in gwp.keys():
	org = lib.split('_')[0]
	if not org in gwpeaks.keys():
		gwpeaks[org] = {}
	for tx in gwp[lib].keys():
		if not tx in gwpeaks[org].keys():
			gwpeaks[org][tx] = []
		gwpeaks[org][tx].append(gwp[lib][tx])

# unique peaks
for org in gwpeaks.keys():
	for tx in gwpeaks[org].keys():
		gwpeaks[org][tx] = list(set(gwpeaks[org][tx]))


## FIND CONSERVED STALL SITES (conserved_stall_sites.py)
organisms = ['yeast', 'fruitfly', 'zebrafish', 'mouse', 'human']

# read BioMart homologs
ens_gname = readBiomart(organisms[1:])
ens_gname['yeast'] = {}
export = open(os.path.join(path, 'biomart_export', 'yeast.txt'), 'r')
for line in export:
	line = line.strip()
	tx = line.split(',')[0]
	gene = line.split(',')[1]
	if not gene == 'n/a':
		if not gene == '':
			if tx in gwpeaks['yeast'].keys():
				if not gene in ens_gname['yeast']:
					ens_gname['yeast'][gene] = [tx]
				else:
					ens_gname['yeast'][gene].append(tx)

# find conserved peaks
conSS = {} ##### !!!!!!!!!!
lostSS = {}

for n in range(len(organisms), 1, -1):
	for curr_org_list in combinations(organisms, n):
		print('gene present in: ', curr_org_list)
		### GET GENES IN COMMON
		gene_tr_list = [ens_gname[org] for org in curr_org_list] # is it in organisms or curr_org_list?
		genes_in_common = set.intersection(*map(set, gene_tr_list)) # set
		for gene in genes_in_common:
			if not gene in conSS: ##### !!!!!!!!!!
				conSS[gene] = {}
				lostSS[gene] = {}
			### get longest transcript
			### write to homologs.fasta
			homologs = open('homologs.fasta', 'w')
			for org in curr_org_list:
				fasta_list = []
				for i in range(len(ens_gname[org][gene])):
					if ens_gname[org][gene][i] in fasta[org]:
						fasta_list.append(fasta[org][ens_gname[org][gene][i]])
					else:
						fasta_list.append('X') ### if fasta_list == ['X'], do not use this!!!
				max_index = fasta_list.index(max(fasta_list, key=len))
				homologs.write('>'+org+'_'+ens_gname[org][gene][max_index]+'\n'+fasta[org][ens_gname[org][gene][max_index]]+'\n')
			homologs.close()
			### run ClustalO
			clustalo = subprocess.Popen(['clustalo', '-i', 'homologs.fasta'], stdout=subprocess.PIPE)
			# to save alignment files, uncomment the following 2 lines (and add a directory to save aln files)
			# out_file_name = gene+'_'+'_'.join(curr_org_list)+'.txt'
			# subprocess.Popen(['clustalo', '-i', 'homologs.fasta', '-o', out_file_name])
			alignment = parse_clustalo(clustalo)
			### get absolute and relative peaks
			abs_peaks = {}
			rel_peaks = {}
			for key in alignment:
				org, tr = key.split('_')[0], key.split('_')[1]
				abs_peaks[org] = []
				rel_peaks[org] = []
				for peak in gwpeaks[org][tr]:
					abs_peaks[org].append(find_absolute_peak_position(alignment[key], peak))
					rel_peaks[org].append((peak, tr))
			### GET STALL SITES IN COMMON
			for n2 in range(len(curr_org_list), 1, -1):
				for curr_org_list2 in combinations(curr_org_list, n2):
					print('stall site conserved in: ', curr_org_list2)
					############# GET CODE FOR PEAKS FROM BELOW HERE !!!!
					abs_peaks_list = [(abs_peaks[org], org) for org in curr_org_list2]
					for peak in abs_peaks_list[0][0]:
						indices = {}
						indices[abs_peaks_list[0][1]] = abs_peaks_list[0][0].index(peak)
						conSS[gene][str(peak)] = []
						conSS[gene][str(peak)].append([abs_peaks_list[0][1], rel_peaks[abs_peaks_list[0][1]][indices[abs_peaks_list[0][1]]]]) ##### !!!!!!!!!!
						lostSS[gene][str(peak)] = []
						for i in range(1, len(abs_peaks_list)):
							ps = abs_peaks_list[i][0]
							if ps:
								if peak-3 <= min(ps, key=lambda x:abs(x-peak)) <= peak+3:
									indices[abs_peaks_list[i][1]] = abs_peaks[abs_peaks_list[i][1]].index(min(ps, key=lambda x:abs(x-peak)))
									conSS[gene][str(peak)].append([abs_peaks_list[i][1], rel_peaks[abs_peaks_list[i][1]][indices[abs_peaks_list[i][1]]]]) ##### !!!!!!!!!!
						if len(indices) == len(abs_peaks_list):
							print('have a match:', peak) # + DEL PEAK
							print([(org, rel_peaks[org][indices[org]][0]) for org in indices])
							for org in abs_peaks:
								if org in indices:
									del abs_peaks[org][indices[org]]
									del rel_peaks[org][indices[org]]
						else:
							del conSS[gene][str(peak)]
					### push the org with missing stall site to separate dictionary
						for o in curr_org_list:
							if not o in curr_org_list2 and len(indices) == len(abs_peaks_list): #### AND IF REMAINING ORGANISMS HAVE THE STALL SITE CONSERVED!
								if rel_peaks[o]:
									k = str(o)+'_'+str(rel_peaks[o][0][1])
									seq = alignment[k][0:peak]
									rel_p = len(seq) - seq.count('-')
									lostSS[gene][str(peak)].append([o, (rel_p, rel_peaks[o][0][1])])
						if not lostSS[gene][str(peak)]:
							del lostSS[gene][str(peak)]
			for org in ens_gname:
				if gene in ens_gname[org]:
					del ens_gname[org][gene] ### DELETE
			if not conSS[gene]:
				del conSS[gene]
			if not lostSS[gene]:
				del lostSS[gene]


# reformat data
conserved_stall_sites = {}
for gene in conSS:
	conserved_stall_sites[gene] = {}
	for peak in conSS[gene]:
		conserved_stall_sites[gene][peak] = {}
		for org in conSS[gene][peak]:
			conserved_stall_sites[gene][peak][org[0]] = org[1]


lost_stall_sites = {}
for gene in lostSS:
	lost_stall_sites[gene] = {}
	for peak in lostSS[gene]:
		lost_stall_sites[gene][peak] = {}
		for org_tr in lostSS[gene][peak]:
			lost_stall_sites[gene][peak][org_tr[0]] = org_tr[1]



### SAVE/EXPORT DATA (save_data.py)

# save_data(conserved_stall_sites, filepath='conserved_stall_sites.p.gz')
# conserved_stall_sites = load_data('conserved_stall_sites.p.gz')

# export to R for plotting and analyses
export_to_R_data_frame(conserved_stall_sites, filename='conserved_stall_sites.csv')


