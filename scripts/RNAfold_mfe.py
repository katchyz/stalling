### fold mRNAs with stall sites
import sys
import RNA

#fin = sys.argv[1]
fout = open(sys.argv[2], 'w')

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
		fasta[header_cds] = seq_cds
	return fasta

fasta = get_FASTA_sequence(sys.argv[1])

for tx in fasta.keys():
	# fout.write('>'+tx+'\n')
	for i in range(0,len(fasta[tx])-50):
		(ss, mfe) = RNA.fold(fasta[tx][i:i+51])
		fout.write(tx+'\t'+str(mfe)+'\n')


fout.close()

