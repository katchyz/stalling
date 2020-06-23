### export stall site conservation (for UpSetR)

import gzip
import cPickle as pickle

ss = pickle.load(gzip.open("conserved_stall_sites.p.gz", "rb"))

organisms = ["yeast", "fruitfly", "zebrafish", "mouse", "human"]

f = open("upsetr.csv", "w")
f.write('gene_ss'+','+','.join([str(x) for x in organisms])+'\n')

for key in ss.keys():
	for key2 in ss[key].keys():
		f.write(key+'_'+key2+',')
		consorg = ss[key][key2].keys()
		#[org in d for org in organisms]
		boo = [int(elem) for elem in [org in consorg for org in organisms]]
		f.write(','.join([str(x) for x in boo]))
		f.write('\n')


f.close()
