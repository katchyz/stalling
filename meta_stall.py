import os
import cPickle as pickle
import gzip
import numpy as np
from scipy import sparse, stats



from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
orfs = BED_get_CDS_mRNA_coord('/Users/kasia/Documents/PhD/data/zebrafish.bed')
(intervals, geneToInterval, geneLength) = BED_to_intervals('/Users/kasia/Documents/PhD/data/zebrafish.bed')



stages = {}
directory = '/Users/kasia/Documents/PhD/meta_data/split_by_length/sparse'
files = filter(lambda f: not f.startswith('.'), os.listdir(directory))
for file in files:
	stages[file[2:-5]] = pickle.load(gzip.open(os.path.join(directory, file), 'rb'))

# One Transcript Per Gene - use only 26459
otpg = {}
f = open('../meta_data/OneTrPerGene.txt', 'r')
for line in f:
	otpg[line.rstrip()] = 1

f.close()

tr_gene = {}
tg = open('/Users/kasia/Documents/PhD/data/zebrafish_tr_gene.txt', 'r')
for line in tg:
	t = line.split()[0]
	g = line.split()[1]
	tr_gene[t] = g

gene_tr = {}
tg = open('/Users/kasia/Documents/PhD/data/zebrafish_tr_gene.txt', 'r')
for line in tg:
	t = line.split()[0]
	g = line.split()[1]
	if g in gene_tr.keys():
		gene_tr[g].append(t)
	else:
		gene_tr[g] = []
		gene_tr[g].append(t)


codon_cov = {}
for stage in stages:
	codon_cov[stage] = {}

codon_cov_high = {}
for stage in stages:
	codon_cov_high[stage] = {}


for stage in stages:
	for tr in otpg: # TR in OTPG
		# merge different read lengths for CDS
		nt_cov = np.asarray(stages[stage][tr][1].sum(axis=0))[0]
		# calculate codon coverage
		if (len(nt_cov) % 3 == 0) and (geneLength[tr] > 33):
			pep_array = np.zeros(len(nt_cov)/3)
			for i in range(0, len(nt_cov), 3):
				pep_array[i/3] = nt_cov[i] + nt_cov[i+1] + nt_cov[i+2]
			codon_cov[stage][tr] = pep_array
			# choose CDSs where at least 50% of codons have at least one footprint
			if np.median(pep_array) > 0:
				codon_cov_high[stage][tr] = pep_array


zscores = {}
for stage in stages:
	zscores[stage] = {}

# get z-scores
## exclude first and last codon in z-score calculation [1:-1]
## two last codons: stop codon and last translating codon
for stage in codon_cov_high:
	for tr in codon_cov_high[stage]:
		zs = stats.zscore(codon_cov_high[stage][tr][1:-2])
		zs = np.insert(zs, 0, 0) # add a zero for 1st codon
		zs = np.append(zs, [0, 0]) # add a zero for last two codons
		zscores[stage][tr] = zs
		# zscores[stage][tr] = stats.zscore(codon_cov_high[stage][tr])

threshold = 8.0

# get genes with peaks (save peak position respective to nucleotide sequence on CDS)
genes_with_peaks = {}
for stage in stages:
	genes_with_peaks[stage] = {}

for stage in zscores:
	for tr in zscores[stage]:
		if max(zscores[stage][tr]) > threshold:
			peak_positions = []
			for i in range(len(zscores[stage][tr])):
				if zscores[stage][tr][i] > threshold:
					peak_positions += [i*3] # 1st nt in codon
			genes_with_peaks[stage][tr] = peak_positions





###########################

meta_stall = {}
for stage in stages:
	meta_stall[stage] = np.zeros((7,33))

for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		for peak in genes_with_peaks[stage][tr]:
			start = peak-15
			end = peak+18
			for l in range(7):
				if start < 0:
					meta_stall[stage][l][abs(start):] += stages[stage][tr][1].toarray()[l][0:end]
				elif end > len(stages[stage][tr][1].toarray()[l]):
					meta_stall[stage][l][0:33-(end-len(stages[stage][tr][1].toarray()[l]))] += stages[stage][tr][1].toarray()[l][start:len(stages[stage][tr][1].toarray()[l])]
				else:
					meta_stall[stage][l] += stages[stage][tr][1].toarray()[l][start:end]

##############

############################
# read in FASTA sequences

fcds = gzip.open('/Users/kasia/Documents/PhD/data/FASTA/cds/Danio_rerio.Zv9.cds.all.fa.gz', 'rb')
header_cds = ''
seq_cds = ''
fasta = {}
for line in fcds:
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


fpep = gzip.open('/Users/kasia/Documents/PhD/data/FASTA/pep/Danio_rerio.Zv9.pep.all.fa.gz', 'rb')
header_pep = ''
seq_pep = ''
peptide = {}
for line in fpep:
	line = line.strip()
	if line.startswith('>'):
		if header_pep and seq_pep:
			peptide[header_pep] = seq_pep
		header_pep = line.split()[4][11:]
		seq_pep = ''
	else:
		seq_pep += line
	# for the last sequence:
	peptide[header_pep] = seq_pep

############

### DNA ###
dna_logo = {}
for stage in stages:
	dna_logo[stage] = []

for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		for peak in genes_with_peaks[stage][tr]:
			start = peak-15
			end = peak+18
			if start < 0:
				dna_logo[stage].append('N'*abs(start)+fasta[tr][0:end])
			elif end > len(fasta[tr]):
				dna_logo[stage].append(fasta[tr][start:len(fasta[tr])]+'N'*(end-len(fasta[tr])))
			else:
				dna_logo[stage].append(fasta[tr][start:end])


f = open('../meta_data/dna_stall-15+18_OTPG.txt', 'w')
for stage in dna_logo:
	f.write('>'+stage+'\n')
	for seqlist in dna_logo[stage]:
		f.write(seqlist+'\n')

### PEPTIDE ###
pep_logo = {}
for stage in stages:
	pep_logo[stage] = []

for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		for peak in genes_with_peaks[stage][tr]:
			start = (peak/3)-50
			end = (peak/3)+11
			if start < 0:
				pep_logo[stage].append('X'*abs(start)+peptide[tr][0:end])
			elif end > len(peptide[tr]):
				pep_logo[stage].append(peptide[tr][start:len(peptide[tr])]+'N'*(end-len(peptide[tr])))
			else:
				pep_logo[stage].append(peptide[tr][start:end])


f = open('/Users/kasia/Documents/PhD/meta_data/peptide_stall_X.txt', 'w')
for stage in pep_logo:
	f.write('>'+stage+'\n')
	for seqlist in pep_logo[stage]:
		f.write(seqlist+'\n')




# check for consistency accross peaks (number of occurences)
# dict of transcripts, each points to two lists: first with position, second with count of occurences at that position
pcons = {}
for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		if tr not in pcons.keys():
			pcons[tr] = [], []
		for i in range(len(genes_with_peaks[stage][tr])):
			if genes_with_peaks[stage][tr][i] in pcons[tr][0]:
				#pcons[tr][0].append(genes_with_peaks_uniq[stage][tr][i])
				pcons[tr][1][pcons[tr][0].index(genes_with_peaks[stage][tr][i])] += 1
			else:
				pcons[tr][0].append(genes_with_peaks[stage][tr][i])
				pcons[tr][1].append(1)
			# counter how many of each...


# count occurences of 1, 2... 8 in pcons[tr][1]
cons = [0]*8
for tr in pcons:
	for i in range(len(pcons[tr][1])):
		cons[pcons[tr][1][i]-1] += 1


# look for repeated peaks a number of times accross stages:
for tr in pcons:
	if 7 in pcons[tr][1]:
		print tr
################
# ten transcripts with peaks in all stages:
# ENSDART00000062324 447	h2afvI	last translating codon
# ENSDART00000008454 489	skp1	last translating codon
# ENSDART00000013575 1254	bzw1a	last translating codon
# ENSDART00000018255 1086	ilf2	double proline
# ENSDART00000011699 192	nono	Pro & others (?)
# ENSDART00000135207 300	rbm4.3  v.high peak, has 2 out of 3 G's downstream
# ENSDART00000056295 1557	psap	last translating codon
# ENSDART00000147057 582	cdh1	has some G's downstream
# ENSDART00000124736 306	HIST1H4F	last translating codon
# ENSDART00000092906 288	ppp1cab		??

# 19 transcripts with peaks in 7 out of 8 stages:
# ENSDART00000125633	khdrbs1a	v.highs, has 4 G's downstream
# ENSDART00000058913	eif4a1a		Glu-Glu, some G's downstream
# ENSDART00000109665	ppm1g		last translating codon
# ENSDART00000111160	ptmab		ON A SCCAFOLD
# ENSDART00000023834	setb 		last translating codon
# ENSDART00000128717	ran 		on a splice site, some G's downstream
# ENSDART00000018475	snrpd3		last translating codon
# ENSDART00000012630	prmt1		last translating codon
# ENSDART00000137689	hnrnpm		last translating codon (not in IGV)
# ENSDART00000002932	marcksb		last translating codon
# ENSDART00000139385	si:ch211-288g17.3 	last translating codon (not in IGV)
# ENSDART00000141703	snrpa		??
# ENSDART00000009790	cx43.4 		has some G's downstream
# ENSDART00000143340	cbx5		last translating codon
# ENSDART00000051799	hdac1		GAA's downstream, out-of-frame stop codon
# ENSDART00000098408	hnrnpabb in IGV: zgc:77052	??
# ENSDART00000112394	ywhabb		2 peaks, plenty of G's downstream
# ENSDART00000135776	srsf7a in IGV: zgc:77155	last translating codon
# ENSDART00000049109	seta		last translating codon

# write peaks to BED file
bed = open('../outputs/stall_sites_OTPG.bed', 'w')
# first line:
# track name="Stall sites" visibility=2 itemRgb="On"
# or:
# colorByStrand="200,100,0 0,100,200"
bed.write('track name="Stall sites" visibility=2 itenRgb="On"\n')
# then:
# chr#	st 	end 	Shield	0	+	(st 	end 	rgb)

for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		# stitch the geneToInterval intervals into mRNA with genomic positions
		mRNA = []
		for exon in range(len(geneToInterval[tr])):
			mRNA.extend(range(geneToInterval[tr][exon].start, geneToInterval[tr][exon].end))
		# offset by CDS start position
		if geneToInterval[tr][0].strand == '+':
			CDS_offset = orfs[tr][0][0]
		elif geneToInterval[tr][0].strand == '-':
			CDS_offset = geneLength[tr] - orfs[tr][0][1] + 3
		for i in range(len(genes_with_peaks[stage][tr])):
			# offset by genes_with_peaks_uniq[stage][tr][i]
			if geneToInterval[tr][0].strand == '+':
				peak_pos = CDS_offset + genes_with_peaks[stage][tr][i]
				# get genomic positions
				gen_st_pos = mRNA[peak_pos]
				gen_end_pos = mRNA[peak_pos]+3 # w/o checking for split over exons...
				# print in bed
				bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'200,100,0'+'\n')
			elif geneToInterval[tr][0].strand == '-':
				peak_pos = CDS_offset + genes_with_peaks[stage][tr][i] - 3
				gen_st_pos = mRNA[-peak_pos]-3
				gen_end_pos = mRNA[-peak_pos]
				bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'0,100,200'+'\n')








#small = orfs[tr][0][0]
#large = orfs[tr][0][1]





##########################
##########################

gwpu_nostend = {}
for stage in stages:
	gwpu_nostend[stage] = {}

for stage in genes_with_peaks:
	for tr in genes_with_peaks[stage]:
		gwpu_nostend[stage][tr] = []
		for i in range(len(genes_with_peaks[stage][tr])):
			if genes_with_peaks[stage][tr][i] != len(stages[stage][tr][1].toarray()[0])-3 and \
				genes_with_peaks[stage][tr][i] != len(stages[stage][tr][1].toarray()[0])-6 and \
				genes_with_peaks[stage][tr][i] != 0:
				gwpu_nostend[stage][tr] += [genes_with_peaks[stage][tr][i]]
		if gwpu_nostend[stage][tr] == []:
			del gwpu_nostend[stage][tr]

for stage in gwpu_nostend:
	print stage
	no = 0
	for tr in gwpu_nostend[stage]:
		no += len(gwpu_nostend[stage][tr])
	print no
# No. peaks in gwpu_nostend
# 28hpf	704
# Shield	1907
# 2-4Cell	110
# 5dpf	528
# 1KCell	265
# Dome 	2598
# 256Cell	200
# Bud 	362

pc_nostend = {}
for stage in gwpu_nostend:
	for tr in gwpu_nostend[stage]:
		if tr not in pc_nostend.keys():
			pc_nostend[tr] = [], []
		for i in range(len(gwpu_nostend[stage][tr])):
			if gwpu_nostend[stage][tr][i] in pc_nostend[tr][0]:
				#pcons[tr][0].append(genes_with_peaks_uniq[stage][tr][i])
				pc_nostend[tr][1][pcons[tr][0].index(gwpu_nostend[stage][tr][i])] += 1
			else:
				pc_nostend[tr][0].append(gwpu_nostend[stage][tr][i])
				pc_nostend[tr][1].append(1)
