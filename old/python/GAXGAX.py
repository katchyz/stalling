# get fasta sequences of (well-expressed) transcripts
# pattern matching for GAXGAXGAX...
# get transcript coverage - check for peak before pattern (z-score? plot distribution of these)


import cPickle as pickle
import gzip
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
from stalling import *
import re

mouse = {}
mouse['mouse_none'] = pickle.load(gzip.open('/Volumes/USELESS/META/all_the_crap/ESC_none_added.p.gz', 'rb'))

bedpath = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

codon_cov_high = get_well_expressed_transcripts(mouse, orfs)

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')

fasta_wet = {}
for tr in codon_cov_high['mouse_none']:
	if tr in fasta:
		fasta_wet[tr] = fasta[tr]


### pattern matching GAXGAXGAXGAXG or GAXGAXXXXXXXG

matches = {}
for tr in fasta_wet:
	match = re.finditer('GA.GA.GA.GA.G', fasta_wet[tr])
	m = []
	for i in match:
		m.append(i.start())
	if m:
		matches[tr] = m

# orf = [4256, 79, 499]
### transcript coverage
zscores = calculate_zscores(codon_cov_high)

del matches['ENSMUST00000129241']
del matches['ENSMUST00000102588']
del matches['ENSMUST00000170213']

z_before = []
for tr in matches:
	print tr
	for pos in matches[tr]:
		print pos
		if pos > 6:
			c = pos/3 - 1
			#z_before.append(max(zscores['mouse_none'][tr][c-1:c+2]))
			print zscores['mouse_none'][tr][c-1:c+2]


### pattern matching  GA.... ... GA.GA.GA.GA.G or .A..A..A.GA.... ... GA.GA.GA.GA.G

matches = {}
for tr in fasta_wet:
	match = re.finditer('.A..A..A.GA..A.GA.GA.GA.GA.GA.G', fasta_wet[tr])
	m = []
	for i in match:
		m.append(i.start())
	if m:
		matches[tr] = m

### transcript coverage


z_psite = []
for tr in matches:
	print tr
	for pos in matches[tr]:
		print pos
		c = pos/3 + 5
		z_psite.append(zscores['mouse_none'][tr][c])
		#print zscores['mouse_none'][tr][c-1:c+2]





