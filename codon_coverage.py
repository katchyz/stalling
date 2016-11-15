### CODON COVERAGE
### mouse, zebrafish

import cPickle as pickle
import gzip
from stalling import get_FASTA_sequence
import numpy as np

# counts
counts_zebrafish = pickle.load(gzip.open("/Volumes/USELESS/META/OUR.p.gz", "rb"))
counts_mouse = pickle.load(gzip.open("/Volumes/USELESS/META/mouse_none.p.gz", "rb"))

# fasta
fasta_zebrafish = get_FASTA_sequence("/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz")
fasta_mouse = get_FASTA_sequence("/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz")

# transcripts per gene
def otpg(gtf_filepath):
	otpg = {}
	gtf = open(gtf_filepath)
	for line in gtf:
		if not line.startswith('#'):
			if line.split()[2] != 'gene': # 1st gene/transcript pair
				gene = line.split()[9][1:-2]
				tr = line.split()[13][1:-2]
				#otpg_mouse[tr] = gene
				if not gene in otpg:
					otpg[gene] = [tr]
				else:
					otpg[gene].append(tr)
	for gene in otpg:
		otpg[gene] = set(otpg[gene])
	return otpg

otpg_zebrafish = otpg('/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf')
otpg_mouse = otpg('/Volumes/USELESS/DATA/genomes/GTF/Mus_musculus.GRCm38.79.chr.gtf')


## get the longest CDS in a set of transcripts

cds_zebrafish = {}
for gene in otpg_zebrafish:
	l = []
	for tr in otpg_zebrafish[gene]:
		if tr in counts_zebrafish['Shield']:
			l.append(counts_zebrafish['Shield'][tr][1].shape[1]) # append length of transcript
		else:
			l.append(0)
	if sum(l) > 0:
		ltr = list(otpg_zebrafish[gene])[np.argmax(l)]
		cds_zebrafish[ltr] = counts_zebrafish['Shield'][ltr][1]

cds_mouse = {}
for gene in otpg_mouse:
	l = []
	for tr in otpg_mouse[gene]:
		if tr in counts_mouse['mouse_none']:
			l.append(counts_mouse['mouse_none'][tr][1].shape[1]) # append length of transcript
		else:
			l.append(0)
	if sum(l) > 0:
		ltr = list(otpg_mouse[gene])[np.argmax(l)]
		cds_mouse[ltr] = counts_mouse['mouse_none'][ltr][1]


## add up all the codons and coverage (per fragment length too)
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]

codons_zebrafish = {}
codons_mouse = {}

for c in codons:
	codons_zebrafish[c] = np.zeros(7)
	codons_mouse[c] = np.zeros(11)



for tr in cds_zebrafish:
	if tr in fasta_zebrafish:
		for i in range(0, cds_zebrafish[tr][1].shape[1], 3):
			cod = fasta_zebrafish[tr][i:i+3]
			if len(cod) == 3 and cod in codons:
				for j in range(7):
					cov = sum(cds_zebrafish[tr][j].todense().tolist()[0][i:i+3])
					codons_zebrafish[cod][j] += cov

# get a number of certain codons too

pickle.dump(codons_zebrafish, gzip.open('/Volumes/USELESS/META/codons_zebrafish.p.gz', 'wb'))

codcount_zebrafish = {}
for c in codons:
	codcount_zebrafish[c] = 0


for tr in cds_zebrafish:
	if tr in fasta_zebrafish:
		for i in range(0, cds_zebrafish[tr][1].shape[1], 3):
			cod = fasta_zebrafish[tr][i:i+3]
			if len(cod) == 3 and cod in codons:
				codcount_zebrafish[cod] += 1


codavg = {}
for c in codsum:
	codavg[c] = codsum[c] / codcount_zebrafish[c]

# dicodon frequencies...




####################################
####################################

### load saved codon sums

import cPickle as pickle
import gzip

plt.ion()

zebra = pickle.load(gzip.open('codons_zebrafish.p.gz', 'rb'))
mouse = pickle.load(gzip.open('codons_mouse.p.gz', 'rb'))

for c in zebra:
	print c, zebra[c][0].sum()

for c in zebra:
	print c, zebra[c][0].sum()/zebra[c][1]


import matplotlib.pyplot as plt
plt.plot(range(25,32), zebra['ATG'][0])


for c in zebra:
	plt.plot(range(25,32), zebra[c][0])



for c in zebra:
	plt.plot(range(25,32), zebra[c][0]/zebra[c][1])

plt.figure()
plt.plot(range(25,32), zebra['GGT'][0]/zebra['GGT'][1])
for c in zebra:
	print c
	raw_input(plt.plot(range(25,32), zebra[c][0]/zebra[c][1]))


for c in mouse:
	plt.plot(range(25,36), mouse[c][0]/mouse[c][1])


plt.figure()
plt.plot(range(25,36), mouse['CGA'][0]/mouse['CGA'][1])
for c in mouse:
	print c
	raw_input(plt.plot(range(25,36), mouse[c][0]/mouse[c][1]))



### mouse3
# stages ESC_36h_HiSeq, ESC_ff_HiSeq, ESC_ff_emet_GA

# see stop codons
plt.figure()
for c in mouse3['ESC_ff_HiSeq']:
	plt.plot(range(25,36), mouse3['ESC_ff_HiSeq'][c][0]/mouse3['ESC_ff_HiSeq'][c][1])

# one codon higher: GAT
plt.figure()
for c in mouse3['ESC_36h_HiSeq']:
	plt.plot(range(25,36), mouse3['ESC_36h_HiSeq'][c][0]/mouse3['ESC_36h_HiSeq'][c][1])

plt.figure()
for c in mouse3['ESC_ff_emet_GA']:
	plt.plot(range(25,36), mouse3['ESC_ff_emet_GA'][c][0]/mouse3['ESC_ff_emet_GA'][c][1])


plt.figure()
plt.plot(range(25,36), mouse3['ESC_36h_HiSeq']['GAT'][0]/mouse3['ESC_36h_HiSeq']['GAT'][1])
for c in mouse3['ESC_36h_HiSeq']:
	print c
	raw_input(plt.plot(range(25,36), mouse3['ESC_36h_HiSeq'][c][0]/mouse3['ESC_36h_HiSeq'][c][1]))

plt.figure()
for c in mouse3['ESC_ff_HiSeq']:
	print c
	raw_input(plt.plot(range(25,36), mouse3['ESC_ff_HiSeq'][c][0]/mouse3['ESC_ff_HiSeq'][c][1]))

plt.figure()
for c in mouse3['ESC_ff_emet_GA']:
	print c
	raw_input(plt.plot(range(25,36), mouse3['ESC_ff_emet_GA'][c][0]/mouse3['ESC_ff_emet_GA'][c][1]))

###################
plt.figure()
plt.plot(range(25,36), mouse['CGA'][0]/mouse['CGA'][1])
for c in mouse:
	print c
	raw_input(plt.plot(range(25,36), mouse[c][0]/mouse[c][1]))


##########################
##########################

## dicodon frequency

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
dicodons = [d+e for d in codons for e in codons]


dicodons_zebrafish = {}
for dic in dicodons:
	dicodons_zebrafish[dic] = 0


for tr in cds_zebrafish:
	print tr
	if tr in fasta_zebrafish:
		for i in range(0, cds_zebrafish[tr][1].shape[1], 3):
			dicod = fasta_zebrafish[tr][i:i+6]
			if len(dicod) == 6 and dicod in dicodons:
				dicodons_zebrafish[dicod] += 1



#####
# get coverage and number of dicodons

for dic in dicodons:
	dicodons_zebrafish[dic] = [np.zeros(7), 0]

for tr in cds_zebrafish:
	if tr in fasta_zebrafish:
		for i in range(0, cds_zebrafish[tr][1].shape[1], 3):
			dicod = fasta_zebrafish[tr][i:i+6]
			if len(cod) == 6 and dicod in dicodons:
				dicodons_zebrafish[dicod][1] += 1
				for j in range(7):
					cov = sum(cds_zebrafish[tr][j].todense().tolist()[0][i:i+6])
					dicodons_zebrafish[dicod][0][j] += cov



import heapq
heapq.nlargest(5, dicodons_zebrafish, key=dicodons_zebrafish.get)












