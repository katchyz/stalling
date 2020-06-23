# translate DNA to amino acids
# make reduced alphabet
import re

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

f = open('yeast.txt', 'r')
o = open('./peptide/yeast.txt', 'w')

for line in f:
	if line.startswith('>'):
		o.write(line)
	else:
		line = line.strip()
		aas = ''
		for i in range(0,len(line),3):
			aas += codon_table[line[i:i+3]]
		o.write(aas)
		o.write('\n')


f.close()
o.close()



f = open('/Volumes/USELESS/STALLING/PLOTS/logo_stall_long/yeast.txt', 'r')
o = open('/Volumes/USELESS/STALLING/PLOTS/logo_stall_long/reduced/yeast.txt', 'w')

positive = re.compile("[RHK]")
negative = re.compile("[DE]")
uncharged = re.compile("[STNQ]")
special = re.compile("[CUGP]")
hydrophobic = re.compile("[AVILMFYW]")
none = re.compile("[X*]")

for line in f:
	if line.startswith('>'):
		o.write(line)
	else:
		line = line.strip()
		aas = ''
		for i in range(0,len(line),3):
			aas += codon_table[line[i:i+3]]
		for aa in aas:
			if positive.match(aa):
				o.write('P')
			elif negative.match(aa):
				o.write('N')
			elif uncharged.match(aa):
				o.write('U')
			elif special.match(aa):
				o.write('S')
			elif hydrophobic.match(aa):
				o.write('H')
			elif none.match(aa):
				o.write('X')
		o.write('\n')




o.close()


