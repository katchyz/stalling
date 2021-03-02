#!/usr/bin/python3

''' translate DNA to amino acids '''

def translate(dna):

	bases = ['t', 'c', 'a', 'g']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))

	protein = ''
	for i in range(0, len(dna), 3):
		if len(dna[i:i+3]) == 3:
			if dna[i:i+3].lower() in codon_table.keys():
				protein += codon_table[dna[i:i+3].lower()]
			else:
				protein += 'N'
		else:
			print('DNA length not multiple of 3, ignoring the overhang')

	return protein


