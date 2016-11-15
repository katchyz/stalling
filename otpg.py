import cPickle as pickle
import gzip
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals, BED_get_CDS_genome_coord
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt

orfs = BED_get_CDS_mRNA_coord('/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed')
(intervals, geneToInterval, geneLength) = BED_to_intervals('/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed')
startpos, endpos = BED_get_CDS_genome_coord('/Volumes/USELESS/DATA/genomes/BED/Danio_rerio.Zv9.79.bed')
# load data (diff. lengths)
our = pickle.load(gzip.open('/Volumes/USELESS/META/OUR.p.gz', 'rb'))

st = '2-4Cell'
# OTPG !!!
otpg = {}
# use one transcript per gene (choose the longest one or random)
# - load GTF
gtf = open('/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf')
f = 0
for line in gtf:
	if not line.startswith('#'):
		if (f == 0) and (line.split()[2] != 'gene'): # 1st gene/transcript pair
			gene = line.split()[9][1:-2]
			tr = line.split()[13][1:-2]
			f = 1
		elif line.split()[2] != 'gene':
			gene_new = line.split()[9][1:-2]
			tr_new = line.split()[13][1:-2]
			if (gene_new == gene) and (tr_new != tr) and (tr in our[st]):
				tr_len = len(our[st][tr][0].toarray()[0])
				tr_new_len = len(our[st][tr_new][0].toarray()[0])
				if tr_new_len > tr_len:
					tr = tr_new
			elif (gene_new != gene):
				### do sth with the existing pair of gene/tr
				otpg[tr] = gene
				gene = gene_new
				tr = tr_new
# dump the last pair
otpg[tr] = gene

# get utr5, cds, utr3
# for different lengths

features = {}

for stage in our:
	features[stage] = {}
	for tr in our[stage]: ### otpg OR our[stage]
		#print tr
		if geneToInterval[tr][0].strand == '+':
			cdsStart = orfs[tr][0][0]
			cdsEnd = orfs[tr][0][1]
		elif geneToInterval[tr][0].strand == '-':
			cdsStart = geneLength[tr] - orfs[tr][0][1]
			cdsEnd = geneLength[tr] - orfs[tr][0][0]
		## check if positive
		if cdsStart >= 0 and cdsEnd >= 0:
			utr5 = our[stage][tr][0][:,0:cdsStart]
			cds = our[stage][tr][0][:,cdsStart:cdsEnd]
			utr3 = our[stage][tr][0][:,cdsEnd::]
			features[stage][tr] = utr5, cds, utr3




pickle.dump(features, gzip.open('/Volumes/USELESS/META/OUR_features.p.gz', 'wb'))
########## for meta_stop...
# if geneToInterval[tr][0].strand == '+':

### plus strand:
# start = orfs[tr][0][1]-18
# end = orfs[tr][0][1]+15
### minus strand:
# start = geneLength[tr]-orfs[tr][0][0]-18
# end = geneLength[tr]-orfs[tr][0][0]+15

# meta_stop[stage][i] += stages[stage][tr][0].toarray()[i][start:end]



# for each length:
# no.fragments per length of feature (store in array)
# (related to the corresponding CDS?)

lct = {}

for stage in features:
	print stage
	lct[stage] = np.zeros([7,3])
	for i in range(7):
		u5, cs, u3 = 0, 0, 0
		for tr in otpg: # otpg OR features[stage]
			if tr in features[stage]:
				u5 += features[stage][tr][0][i].sum()
				cs += features[stage][tr][1][i].sum()
				u3 += features[stage][tr][2][i].sum()
		print "5'UTR: ", u5, "CDS: ", cs, "3'UTR: ", u3
		lct[stage][i] = np.array([u5, cs, u3])



p = 0
plt.figure('leaders')
for stage in features:
	print stage
	for fl in range(0,7):
		cum_ribocov = np.array([])
		for tr in features[stage]:
			if tr in features[stage]:
				ribocov = np.array(features[stage][tr][0].todense()[fl])[0] # 3: len28
				if not cum_ribocov.any(): # 1st round
					cum_ribocov = ribocov
					l = len(ribocov) # ??? only start ???
				elif len(ribocov) > len(cum_ribocov):
					cum_ribocov = np.append(cum_ribocov, np.zeros(len(ribocov)-len(cum_ribocov)))
					cum_ribocov += ribocov
				else: # len(cum_ribocov) >= len(ribocov)
					cum_ribocov[0:len(ribocov)] += ribocov
					l = len(ribocov) # ??? only start ???
		ribocov = cum_ribocov[0:500] # ??? only start ??? # [offset:l] # offset=12 ### GET FIRST 500 NT ### if we have such long CDSs
		nts = np.array(range(0, len(ribocov))) # range(0, len(ribocov)-offset)
		ft_ribocov = np.fft.fft(ribocov, axis=0)
		amplitudes = abs(ft_ribocov)
		frequencies = np.fft.fftfreq(len(ribocov), 1)
		periods = 1 / frequencies
		plt.subplot(10,7,p+fl+1)
		plt.plot(periods, amplitudes * 1e-3)
		plt.xlim(0, 10)
	p += 7


plt.show()


p = 0
plt.figure('trailers')
for stage in features:
	print stage
	for fl in range(0,7):
		cum_ribocov = np.array([])
		for tr in features[stage]:
			if tr in features[stage]:
				ribocov = np.array(features[stage][tr][2].todense()[fl])[0] # 3: len28
				if not cum_ribocov.any(): # 1st round
					cum_ribocov = ribocov
					l = len(ribocov) # ??? only start ???
				elif len(ribocov) > len(cum_ribocov):
					cum_ribocov = np.append(cum_ribocov, np.zeros(len(ribocov)-len(cum_ribocov)))
					cum_ribocov += ribocov
				else: # len(cum_ribocov) >= len(ribocov)
					cum_ribocov[0:len(ribocov)] += ribocov
					l = len(ribocov) # ??? only start ???
		ribocov = cum_ribocov[0:500] # ??? only start ??? # [offset:l] # offset=12 ### GET FIRST 500 NT ### if we have such long CDSs
		nts = np.array(range(0, len(ribocov))) # range(0, len(ribocov)-offset)
		ft_ribocov = np.fft.fft(ribocov, axis=0)
		amplitudes = abs(ft_ribocov)
		frequencies = np.fft.fftfreq(len(ribocov), 1)
		periods = 1 / frequencies
		plt.subplot(10,7,p+fl+1)
		plt.plot(periods, amplitudes * 1e-3)
		plt.xlim(0, 10)
	p += 7


plt.show()


p = 0
plt.figure('CDSs')
for stage in features:
	print stage
	for fl in range(0,7):
		cum_ribocov = np.array([])
		for tr in features[stage]:
			if tr in features[stage]:
				ribocov = np.array(features[stage][tr][1].todense()[fl])[0] # 3: len28
				if not cum_ribocov.any(): # 1st round
					cum_ribocov = ribocov
					l = len(ribocov) # ??? only start ???
				elif len(ribocov) > len(cum_ribocov):
					cum_ribocov = np.append(cum_ribocov, np.zeros(len(ribocov)-len(cum_ribocov)))
					cum_ribocov += ribocov
				else: # len(cum_ribocov) >= len(ribocov)
					cum_ribocov[0:len(ribocov)] += ribocov
					l = len(ribocov) # ??? only start ???
		ribocov = cum_ribocov[0:500] # ??? only start ??? # [offset:l] # offset=12 ### GET FIRST 500 NT ### if we have such long CDSs
		nts = np.array(range(0, len(ribocov))) # range(0, len(ribocov)-offset)
		ft_ribocov = np.fft.fft(ribocov, axis=0)
		amplitudes = abs(ft_ribocov)
		frequencies = np.fft.fftfreq(len(ribocov), 1)
		periods = 1 / frequencies
		plt.subplot(10,7,p+fl+1)
		plt.plot(periods, amplitudes * 1e-3)
		plt.xlim(0, 10)
	p += 7


plt.show()

