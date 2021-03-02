#!/usr/bin/python3

''' Reads FASTA sequence files. '''

import gzip


def get_FASTA_sequence(filepath):
	""" Reads gzipped fasta file """
 
	if filepath.endswith('.gz'):
		f = gzip.open(filepath, 'rt')
	else:
		f = open(filepath, 'rt')
 
	fasta = {}
	header_cds = ''
	seq_cds = ''
 
	c = 0
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if header_cds and seq_cds:
				fasta[header_cds] = {}
				fasta[header_cds]['seq'] = seq_cds
				fasta[header_cds]['gene_id'] = gene_id
				if c == 0:
					fasta[header_cds]['tx_id'] = tx_id
			info = line.split() # header line
			header_cds = info[0][1:]
			seq_cds = ''
			gene_id = info[3][5:]
			if c == 0 and info[4].startswith('transcript'):
				tx_id = info[4][11:]
			else:
				c = 1
		else:
			seq_cds += str(line)
		# for the last sequence:
		fasta[header_cds] = {}
		fasta[header_cds]['seq'] = seq_cds
		fasta[header_cds]['gene_id'] = gene_id
		if c == 0:
			fasta[header_cds]['tx_id'] = tx_id
 
	f.close()
 
	return fasta


######### EXAMPLE

# fasta_cdna = "../DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
# fasta_pep = "../DATA/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz"

# fp = get_FASTA_sequence(fasta_pep)


