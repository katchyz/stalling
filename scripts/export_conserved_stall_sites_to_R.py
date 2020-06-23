### save stall sites with positions on transcripts
import gnureadline
import cPickle as pickle
import gzip
import pandas as pd

f = pickle.load(gzip.open("conserved_stall_sites.p.gz", 'rb'))

organisms = ['yeast', 'fruitfly', 'zebrafish', 'mouse', 'human']

# df[gene] = [yeast_tx, yeast_pos, fruitfly_tx, fruitfly_pos, zebrafish_tx, zebrafish_pos, mouse_tx, mouse_pos, human_tx, human_pos]
df = {}


for gene in f.keys():
	for ss in f[gene].keys():
		df[gene+'_'+ss] = ['NA']*10 # empty array
		for org in f[gene][ss].keys():
			pos = f[gene][ss][org][0]
			tx = f[gene][ss][org][1]
			# add to df
			df[gene+'_'+ss][organisms.index(org)*2] = tx
			df[gene+'_'+ss][organisms.index(org)*2+1] = pos


pdf = pd.DataFrame.from_dict(df, orient='index')
pdf.to_csv('conserved_stall_sites.csv', header=False)

### add names of (longest) transcripts if homologs are present in other organisms
homologs = {}

for n in range(len(organisms), 1, -1):
	for curr_org_list in combinations(organisms, n):
		print 'gene present in: ', curr_org_list
		### GET GENES IN COMMON
		gene_tr_list = [ens_gname[org] for org in curr_org_list] # is it in organisms or curr_org_list?
		genes_in_common = set.intersection(*map(set, gene_tr_list)) # set
		for gene in genes_in_common:
			if not gene in homologs:
				homologs[gene] = ['NA']*5
			### get longest transcript
			### write to homologs.fasta
			for org in curr_org_list:
				fasta_list = []
				for i in range(len(ens_gname[org][gene])):
					if ens_gname[org][gene][i] in fasta[org]:
						fasta_list.append(fasta[org][ens_gname[org][gene][i]])
					else:
						fasta_list.append('X') ### if fasta_list == ['X'], do not use this!!!
				max_index = fasta_list.index(max(fasta_list, key=len))
				homologs[gene][organisms.index(org)] = ens_gname[org][gene][max_index]

pickle.dump(homologs, gzip.open('homologs.p.gz', 'wb'))

###############
df = pickle.load(gzip.open('homologs.p.gz', 'rb'))
pdf = pd.DataFrame.from_dict(df, orient='index')
pdf.to_csv('homologs.csv', header=False)


