### check for complementarity between 18S rRNA and mRNA fragments around stall sites
### try 28S too
import gzip
import cPickle as pickle
from stalling import get_sequence_logo_around_stall_sites
import random
from shoelaces_full_distribution import BED_to_intervals
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import

# get 18S rRNA sequence

def get_FASTA_sequence(filepath):
	""" Reads gzipped fasta file """
	f = open(filepath, 'r')
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
		fasta[header_cds] = seq_cds
	return fasta

rRNA = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/rRNA/18S.fasta')
rRNA = rRNA['gi|374088232|ref|NR_003278.3|']

# get fragments around stall sites (varying lengths)
css = pickle.load(gzip.open('/Volumes/USELESS/META/conserved_stall_sites.p.gz', 'rb'))

seqs = {}
seqs['mouse'] = {}
for gene in css:
	for peak in css[gene]:
		if 'mouse' in css[gene][peak]:
			seqs['mouse'][css[gene][peak]['mouse'][1]] = [css[gene][peak]['mouse'][0]]

# get the sequence on the transcript...
from stalling import get_FASTA_sequence
fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')

SSseq = get_sequence_logo_around_stall_sites(fasta, seqs, -15, 15)

# (reverse) complement - check sense and antisense (?)
def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
	return reverse_complement


for i in SSseq['mouse']:
	print '>'+str(random.random()), '\n', reverse_complement(i)
# allow mismatches? wobble base pairing? Bulging (there's not much space)



# get control - random fragements from genes
contr = {}
contr['mouse'] = {}
for tr in seqs['mouse']:
	contr['mouse'][tr] = [seqs['mouse'][tr][0] - 50]

CONTRseq = get_sequence_logo_around_stall_sites(fasta, contr, -15, 15)


for i in CONTRseq['mouse']:
	if len(i) == 33:
		print '>'+str(random.random()), '\n', reverse_complement(i)

##########

contr2 = {}
contr2['mouse'] = {}
for tr in seqs['mouse']:
	contr2['mouse'][tr] = [random.randrange(20, len(fasta[tr])-20)]

CONTRseq2 = get_sequence_logo_around_stall_sites(fasta, contr2, -15, 15)

for i in CONTRseq2['mouse']:
	if len(i) == 33:
		print '>'+str(random.random()), '\n', reverse_complement(i)



#############
#############
kmers = []
for s in SSseq['mouse']:
	seq = s[15:28]
	#print seq
	#print '>'+str(random.random()), '\n', reverse_complement(seq)
	kmers.append(seq)


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


pos_18S = [0]*len(rRNA)
for k in kmers:
	krc = reverse_complement(k)
	for i in range(0, len(rRNA)-13):
		seed = rRNA[i:i+13]
		if similar(krc, seed) > 0.7:
			pos_18S[i] += 1

plt.plot(range(len(rRNA)), pos_18S, 'o')
plt.show()

###
pos_18S = [0]*len(rRNA)
for k in kmers:
	krc = reverse_complement(k)
	for i in range(0, len(rRNA)-13):
		seed = rRNA[i:i+13]
		pos_18S[i] += similar(krc, seed)


plt.plot(range(len(rRNA)), pos_18S, 'o')
plt.show()



for k in kmers:
	print '>'+str(random.random()), '\n', reverse_complement(k)

### score seeds from 18S against hmm.out profile (on Desktop)



GA_content = [0]*len(kmers)
for i in range(len(kmers)):
	GA_content[i] += kmers[i].count('G') + kmers[i].count('A')

CT_content = [0]*(len(rRNA)-13)
for i in range(0, len(rRNA)-13):
	seed = rRNA[i:i+13]
	CT_content[i] += seed.count('C') + seed.count('T')





controls2 = []
for seq in CONTRseq2['mouse']:
	controls2.append(seq[18:31])

GA_contr = [0]*len(controls2)
for i in range(len(controls2)):
	GA_contr[i] += controls2[i].count('G') + controls2[i].count('A')



controls = []
for seq in CONTRseq['mouse']:
	controls.append(seq[18:31])

GA_contr1 = [0]*len(controls)
for i in range(len(controls)):
	GA_contr1[i] += controls[i].count('G') + controls[i].count('A')



# try on bigger datasets, more controls (make sure the controls are out of stall site regions)
#############
##############

fragments = []
for s in SSseq['mouse']:
	fragments.append(s)

GA_content = [0]*len(fragments)
for i in range(len(fragments)):
	GA_content[i] += fragments[i].count('G') + fragments[i].count('A')

plt.hist(GA_content, histtype='step')

### many times

contr2 = {}
contr2['mouse'] = {}
for tr in seqs['mouse']:
	contr2['mouse'][tr] = [random.randrange(20, len(fasta[tr])-20)]

CONTRseq2 = get_sequence_logo_around_stall_sites(fasta, contr2, -15, 15)

controls2 = []
for seq in CONTRseq2['mouse']:
	#controls2.append(seq[18:31])
	controls2.append(seq)

GA_contr = [0]*len(controls2)
for i in range(len(controls2)):
	GA_contr[i] += controls2[i].count('G') + controls2[i].count('A')

plt.hist(GA_contr, histtype='step')
