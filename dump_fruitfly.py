import numpy as np
import cPickle as pickle
import gzip

def csr_vappend(a,b):
	""" Takes in 2 csr_matrices and appends the second one to the bottom of the first one. 
	Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
	the first matrix instead of copying it. The data, indices, and indptr still get copied."""
	a.data = np.hstack((a.data,b.data))
	a.indices = np.hstack((a.indices,b.indices))
	a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
	a._shape = (a.shape[0]+b.shape[0],b.shape[1])
	return a


d = {}

# stage_dir_names = ['0-2h_embryo_A', '0-2h_embryo_B', '0-2h_embryo_cushion_A', '0-2h_embryo_cushion_B', \
# '0-2h_embryo_gradient_A', '0-2h_embryo_gradient_B', 'S2Cell_150_A', 'S2cell_150_B', 'S2cell_250_A', 'S2cell_250_B']

stage_dir_28_37 = ['0-2h_embryo_A', '0-2h_embryo_B']
stage_dir_28_36 = ['0-2h_embryo_cushion_A', '0-2h_embryo_cushion_B', '0-2h_embryo_gradient_A', '0-2h_embryo_gradient_B']
stage_dir_27_37 = ['S2cell_150_A', 'S2cell_150_B']
stage_dir_27_35 = ['S2cell_250_A', 'S2cell_250_B']

for stage in stage_dir_28_37:
	d[stage] = {}
	for l in range(28, 38):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/fruitfly_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

for stage in stage_dir_28_36:
	d[stage] = {}
	for l in range(28, 37):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/fruitfly_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

for stage in stage_dir_27_37:
	d[stage] = {}
	for l in range(27, 38):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/fruitfly_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

for stage in stage_dir_27_35:
	d[stage] = {}
	for l in range(27, 36):
		d[stage][str(l)] = pickle.load(gzip.open('/Volumes/USELESS/META/fruitfly_meta/'+stage+'-'+str(l)+'.p.gz', 'rb'))

#

seq_orf = {}

##
for stage in stage_dir_28_37:
	seq_orf[stage] = {}
	for tr in d[stage]['28']:
		seq_orf[stage][tr] = d[stage]['28'][tr][0], d[stage]['28'][tr][1]

for stage in stage_dir_28_37:
	for l in range(29, 38):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])

##
for stage in stage_dir_28_36:
	seq_orf[stage] = {}
	for tr in d[stage]['28']:
		seq_orf[stage][tr] = d[stage]['28'][tr][0], d[stage]['28'][tr][1]

for stage in stage_dir_28_36:
	for l in range(29, 37):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])

##
for stage in stage_dir_27_37:
	seq_orf[stage] = {}
	for tr in d[stage]['27']:
		seq_orf[stage][tr] = d[stage]['27'][tr][0], d[stage]['27'][tr][1]

for stage in stage_dir_27_37:
	for l in range(28, 38):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])

##
for stage in stage_dir_27_35:
	seq_orf[stage] = {}
	for tr in d[stage]['27']:
		seq_orf[stage][tr] = d[stage]['27'][tr][0], d[stage]['27'][tr][1]

for stage in stage_dir_27_35:
	for l in range(28, 36):
		for tr in d[stage][str(l)]:
			csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
			csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])


pickle.dump(seq_orf, gzip.open('/Volumes/USELESS/META/fruitfly.p.gz', 'wb'))