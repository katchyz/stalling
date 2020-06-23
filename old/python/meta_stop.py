import os
import cPickle as pickle
import gzip
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
import numpy as np


# get all dictionaries for all stages
stages = {}
directory = '/Users/kasia/Documents/PhD/meta_data/split_by_length/sparse'
files = filter(lambda f: not f.startswith('.'), os.listdir(directory))
for file in files:
	stages[file[2:-5]] = pickle.load(gzip.open(os.path.join(directory, file), 'rb'))

# get CDS coordinates on mRNA
orfs = BED_get_CDS_mRNA_coord('/Users/kasia/Documents/PhD/data/zebrafish.bed')
(intervals, geneToInterval, geneLength) = BED_to_intervals('/Users/kasia/Documents/PhD/data/zebrafish.bed')

# dictionary for eight different stages and array for seven different lengths
meta_stop = {}
for stage in stages:
	meta_stop[stage] = np.zeros((7,33))


# otpg = {}
# f = open('../meta_data/OneTrPerGene.txt', 'r')
# for line in f:
# 	otpg[line.rstrip()] = 1

for tr in orfs: ## OTPG !!!!!!!
	# check if lenght of coding sequence is a multiple of 3 nt
	if (len(stages['2-4Cell'][tr][1].toarray()[0]) % 3 == 0) and (geneLength[tr] > 33):
		if geneToInterval[tr][0].strand == '+':
			for stage in stages:
				for i in range(7):
					start = orfs[tr][0][1]-18
					end = orfs[tr][0][1]+15
					if start < 0:
						meta_stop[stage][i][abs(start):] += stages[stage][tr][0].toarray()[i][0:end]
					elif end > len(stages[stage][tr][0].toarray()[i]):
						meta_stop[stage][i][0:33-(end-len(stages[stage][tr][0].toarray()[i]))] += stages[stage][tr][0].toarray()[i][start:len(stages[stage][tr][0].toarray()[i])]
					else:
						meta_stop[stage][i] += stages[stage][tr][0].toarray()[i][start:end]
		elif geneToInterval[tr][0].strand == '-':
			for stage in stages:
				for i in range(7):
					start = geneLength[tr]-orfs[tr][0][0]-18
					end = geneLength[tr]-orfs[tr][0][0]+15
					if start < 0:
						meta_stop[stage][i][abs(start):] += stages[stage][tr][0].toarray()[i][0:end]
					elif end > len(stages[stage][tr][0].toarray()[i]):
						meta_stop[stage][i][0:33-(end-len(stages[stage][tr][0].toarray()[i]))] += stages[stage][tr][0].toarray()[i][start:len(stages[stage][tr][0].toarray()[i])]
					else:
						meta_stop[stage][i] += stages[stage][tr][0].toarray()[i][start:end]




#small = orfs[tr][0][0]
#large = orfs[tr][0][1]

fcdna = gzip.open('/Users/kasia/Documents/PhD/data/FASTA/cdna/Danio_rerio.Zv9.cdna.all.fa.gz', 'rb')
header_cdna = ''
seq_cdna = ''
cdna = {}
for line in fcdna:
	line = line.strip()
	if line.startswith('>'):
		if header_cdna and seq_cdna:
			cdna[header_cdna] = seq_cdna
		header_cdna = line.split()[0][1:]
		seq_cdna = ''
	else:
		seq_cdna += line
	# for the last sequence:
	cdna[header_cdna] = seq_cdna


cdna_logo = {}
for stage in stages:
	cdna_logo[stage] = []


for tr in otpg:
	# check if lenght of coding sequence is a multiple of 3 nt
	if (len(stages['2-4Cell'][tr][1].toarray()[0]) % 3 == 0) and (geneLength[tr] > 33):
		if geneToInterval[tr][0].strand == '+':
			for stage in stages:
				start = orfs[tr][0][1]-18
				end = orfs[tr][0][1]+15
				if start < 0:
					cdna_logo[stage].append('N'*abs(start)+cdna[tr][0:end])
				elif end > len(cdna[tr]):
					cdna_logo[stage].append(cdna[tr][start:len(cdna[tr])]+'N'*(end-len(cdna[tr])))
				else:
					cdna_logo[stage].append(cdna[tr][start:end])
		elif geneToInterval[tr][0].strand == '-':
			for stage in stages:
				start = geneLength[tr]-orfs[tr][0][0]-18
				end = geneLength[tr]-orfs[tr][0][0]+15
				if start < 0:
					cdna_logo[stage].append('N'*abs(start)+cdna[tr][0:end])
				elif end > len(cdna[tr]):
					cdna_logo[stage].append(cdna[tr][start:len(cdna[tr])]+'N'*(end-len(cdna[tr])))
				else:
					cdna_logo[stage].append(cdna[tr][start:end])


f = open('../meta_data/cdna_stop-18+15_OTPG.txt', 'w')
for stage in cdna_logo:
	f.write('>'+stage+'\n')
	for seqlist in cdna_logo[stage]:
		f.write(seqlist+'\n')




