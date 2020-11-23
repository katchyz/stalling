#!/usr/bin/python3

''' Finds the longest transcript per gene (optionally),
	selects transcripts with good coverage (median codon coverage > 0),
	calculates z-score coverage,
	calls peaks for each transcript.
'''

import numpy as np
from scipy import stats



def get_longest_transcript(transcripts):
	""" Takes output from read_BED or read_GTF and returns only the longest tx (cds) per gene """
 
	longest_tx = {}
 
	for tx in transcripts:
		gene = transcripts[tx]['gene_id']
		cds_len = transcripts[tx]['cds_coord'][1] - transcripts[tx]['cds_coord'][0]
 
		if not gene in longest_tx.keys():
			longest_tx[gene] = transcripts[tx]
			longest_tx[gene]['tx_id'] = tx
		else:
			if cds_len > longest_tx[gene]['cds_coord'][1] - longest_tx[gene]['cds_coord'][0]:
				longest_tx[gene] = transcripts[tx]
				longest_tx[gene]['tx_id'] = tx
 
	return longest_tx


def set_keys_to_tx_id(longest_tx):
	""" Changes the keys of longest_tx dictionary from gene_id to tx_id. """
 
	longest_tx_id = {}
 
	for gene in longest_tx:
		tx = longest_tx[gene]['tx_id']
		longest_tx_id[tx] = longest_tx[gene]
 
	return longest_tx_id


def get_well_expressed_transcripts(ribo_cov):
	""" Takes a dictionary of ribosome coverage and returns a dictionary with well expressed transcripts
		(median codon coverage less than one). """
 
	well_expr_tx = {}
 
	for tx in ribo_cov:
		# nt_cov = np.asarray(data_dict[stage][tr][1].sum(axis=0))[0] per length
		nt_cov = np.nan_to_num(ribo_cov[tx])
		# calculate codon coverage
		if len(nt_cov) > 33: # len(nt_cov) % 3 == 0
			pep_array = np.zeros(int(np.floor(len(nt_cov)/3)))
			for i in range(0, int(np.floor(len(nt_cov)/3)*3), 3):
				pep_array[int(i/3)] = nt_cov[i] + nt_cov[i+1] + nt_cov[i+2]
 
			# choose CDSs where at least 50% of codons have at least one footprint
			if np.median(pep_array) > 0:
				well_expr_tx[tx] = pep_array
 
		else:
			print(f'Length of gene {tx} is less than 33 nt')
 
	return well_expr_tx


def calculate_zscores(well_expr_tx, start=1, end=2):
	""" Given the dictionary of highly expressed transcripts, calculates
		the z-scores, avoiding first and last codons.
		start and end are numbers of codons at the beginning/end of tx to exclude from calculation. """
 
	zscores = {}
 
	for tx in well_expr_tx:
		zs = stats.zscore(well_expr_tx[tx][start:-end])
		zs = np.insert(zs, 0, [0]*start) # add a zero for 1st codon
		zs = np.append(zs, [0]*end) # add a zero for last two codons
		zscores[tx] = zs
 
	return zscores


def get_genes_with_peaks(zscores, threshold=5):
	""" Given the dictionary with z-scores and a z-score threshold,
		return the dictionary of transcripts with list of peaks. """
 
	genes_with_peaks = {}
 
	for tx in zscores:
		if max(zscores[tx]) > threshold:
			peak_positions = np.where(zscores[tx] > threshold)[0] *3 # 1st nt in codon
			genes_with_peaks[tx] = peak_positions
 
	return genes_with_peaks





longest_tx = get_longest_transcript(transcripts)
longest_tx = set_keys_to_tx_id(longest_tx)
ribo_cov = extract_intervals_from_wigs(wig_fwd, wig_rev, longest_tx)
wet = get_well_expressed_transcripts(ribo_cov)
zscores = calculate_zscores(wet, start=5, end=2)
gwp = get_genes_with_peaks(zscores, threshold=5)





