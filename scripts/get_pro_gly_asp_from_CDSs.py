#!/usr/bin/env python

# get all pro/gly/asp codons in the CDSs, save logos to plot context

import sys

#fin = sys.argv[1]

pro = open('all_pro.txt', 'w')
gly = open('all_gly.txt', 'w')
asp = open('all_asp.txt', 'w')

def get_FASTA_sequence(file):
	""" Reads fasta file """
	f = open(file, 'r')
	fasta = {}
	header_cds = ''
	seq_cds = ''
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if header_cds and seq_cds:
				fasta[header_cds] = seq_cds
			header_cds = line.split()[0][1:]
			seq_cds = ''
		else:
			seq_cds += line
		# for the last sequence:
		seq_cds = seq_cds.upper()
		fasta[header_cds] = seq_cds
	return fasta

fasta = get_FASTA_sequence(sys.argv[1])

for tx in fasta.keys():
	# fout.write('>'+tx+'\n')
	for i in range(15,len(fasta[tx])-20):
		if fasta[tx][i:(i+2)] == "CC":
			seq = fasta[tx][(i-15):(i+18)]
			pro.write(seq+"\n")
		elif fasta[tx][i:(i+2)] == "GG":
			seq = fasta[tx][(i-15):(i+18)]
			gly.write(seq+"\n")
		if fasta[tx][i:(i+3)] == "GAT" or fasta[tx][i:(i+3)] == "GAC":
			seq = fasta[tx][(i-15):(i+18)]
			asp.write(seq+"\n")



pro.close()
gly.close()
asp.close()






