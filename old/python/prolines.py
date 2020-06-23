import re
from stalling import get_FASTA_pep_sequence

fasta_pep = get_FASTA_pep_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/pep/Danio_rerio.Zv9.pep.all.fa.gz')

iterators = {}
for tr in fasta:
	iterators[tr] = re.finditer('P{2,}', fasta_pep[tr])

# n = 0
# for tr in iterators:
# 	for i in list(iterators[tr]):
# 		print tr, i.group(), i.start(), i.end()
# 		n += 1


ps = {}
for tr in iterators:
	ps[tr] = []
	for i in list(iterators[tr]):
		ps[tr].append([i.start(), i.end()])

pp = {}
for tr in ps:
	if not ps[tr] == []:
		pp[tr] = ps[tr]

#########

z = []
for tr in pp:
	if tr in zscores_Giraldez['RPF_05h_column_rep2_run3']:
		for i in pp[tr]:
			# i[0], i[1] - start and end of prolines
			z.append(max(zscores_Giraldez['RPF_05h_column_rep2_run3'][tr][i[0]:i[1]]))

#########

import matplotlib.pyplot as plt
plt.hist(z, len(z))
plt.show()

#########


p = re.finditer('P{2,}', sequence)

for i in list(p):
	print i.group(), i.start(), i.end()

########
########

# look for z-score over any of the P codons (zscores: length of peptide)
# take only the highest one
# plot...

# check how many of these are in our/Giraldez data