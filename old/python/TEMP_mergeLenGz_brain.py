import numpy as np
import cPickle as pickle
import gzip

d = {}
stage_dir_names = ['Brain1_pm', 'Brain2_mn', 'Brain3_am', 'Brain4_nn']

for stage in stage_dir_names:
	d[stage] = {}
	for l in range(25, 32):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/zebrafish_brain_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

def csr_vappend(a,b):
	""" Takes in 2 csr_matrices and appends the second one to the bottom of the first one. 
	Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
	the first matrix instead of copying it. The data, indices, and indptr still get copied."""
	a.data = np.hstack((a.data,b.data))
	a.indices = np.hstack((a.indices,b.indices))
	a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
	a._shape = (a.shape[0]+b.shape[0],b.shape[1])
	return a


#

seq_orf = {}
for stage in d:
	seq_orf[stage] = {}
	for tr in d[stage]['25']:
		seq_orf[stage][tr] = d[stage]['25'][tr][0], d[stage]['25'][tr][1]

for stage in seq_orf:
	for l in range(26, 32):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])


pickle.dump(seq_orf, gzip.open('/Volumes/USELESS/META/Brain.p.gz', 'wb'))