### mRNA structure

# get fasta sequence around stall sites
# otpg ?
# from -100 to +100
# window size of ~50nt?  ### VARY WINDOW SIZE !!

# ... on each sequence ...
# for each window, moving the beginning every 5nt (1nt?) until ~20nt after the stall site
# run RNAfold, store min. free energies

# add them up, make global plot

import RNA
import cPickle as pickle
import gzip
from stalling import *

gwp = pickle.load(gzip.open('/Volumes/USELESS/META/consensus.p.gz', 'rb'))

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz')

def get_sequence(fasta, genes_with_peaks, upstream, downstream):
	""" Gets fasta sequences and list of stall sites, and number of nucleotides or amino acids
	upstream and downstream of the stall codon, returns lists of sequences for each stage.
	Upstream should be a negative integer, downstream a positive """
	seq_logo = {}
	for stage in genes_with_peaks:
		seq_logo[stage] = []
		for tr in genes_with_peaks[stage]:
			if tr in fasta:
				for peak in genes_with_peaks[stage][tr]:
					start = peak + upstream ### vary the length
					end = peak + 3 + downstream ### vary the length
					if start >= 0 and end <= len(fasta[tr]):
						seq_logo[stage].append(fasta[tr][start:end])
	return seq_logo


seq_logo = get_sequence(fasta, gwp, -15, 150)

#########
#########

win_size = 30
energies = []
for seq in seq_logo['consensus']:
	parE = []
	for i in range(0, len(seg_logo['consensus'][0])-win_size, 3):
		win = seq[i:i+win_size]
		parE.append(RNA.fold(win)[1])
	energies.append(parE)



########
import matplotlib.pyplot as plt
x = range(-99, 99-win_size+3, 3)
for y in energies:
	plt.plot(x, y)

plt.show()
########

RNA.fold('CCCCCAAAGCTTTGGGGG')

minE = RNA.fold('CCCATGGG')[1] # float

##########

# check for structure, varying window lengths (just after stall site) - if the structure starts with '('
	# take the lowest minE
# starting some nucleotides before the stall site

# keep the same max range (go only through say ~50 nt, with a window of max 150)

l = len(seq_logo['consensus'][0])
energies = []
for seq in seq_logo['consensus']:
	parE = []
	for i in range(0, l):
		winE = []
		for j in range(5, l-i):
			win = seq[i:i+j]
			s = RNA.fold(win)
			if s[0][0] == '(':
				winE.append(s[1])
		if winE:
			parE.append(min(winE))
		else:
			parE.append('none')
	energies.append(parE)


arrays = []
for i in energies:
	n = np.array([])
	for e in i:
		if type(e) == float:
			n = np.append(n, [e])
		elif e == 'none':
			n = np.append(n, [np.nan])
	arrays.append(n)


import matplotlib.pyplot as plt
Fig, ax = plt.subplots()

x = range(-15, 153)
for y in arrays:
	ax.plot(x,y)

Fig.show()


####################
####################
# keep the same max range (go only through say ~50 nt, with a window of max 150)

l = len(seq_logo['consensus'][0])


energies = []
for seq in seq_logo['consensus']:
	parE = np.array([])
	for i in range(0, 48):
		winE = []
		for j in range(5, 120):
			win = seq[i:i+j]
			s = RNA.fold(win)
			if s[0][0] == '(':
				winE.append(s[1])
		if winE:
			parE = np.append(parE, [min(winE)])
		else:
			parE = np.append(parE, [np.nan])
	energies.append(parE)

pickle.dump(energies, gzip.open('/Volumes/USELESS/META/energies.p.gz', 'wb'))


Fig, ax = plt.subplots()

x = range(-15, 33)
for y in energies:
	ax.plot(x,y)

Fig.show()



