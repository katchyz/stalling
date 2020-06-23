# extract ALL transcripts from wiggle files
# for 8 stages, 7 lengths(25-31)
# get stall positions

from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord, BED_get_CDS_genome_coord


(intervals, geneToInterval, geneLength) = BED_to_intervals('/Users/kasia/Documents/PhD/data/zebrafish.bed')
#set, dict, dict
orfs = BED_get_CDS_mRNA_coord('/Users/kasia/Documents/PhD/data/zebrafish.bed')
#dict
#(startpos, endpos) = BED_get_CDS_genome_coord('/Users/kasia/Documents/PhD/data/zebrafish.bed')
#set, set


import os
from wiggle import ParseWig
import numpy as np
from scipy import sparse
import cPickle as pickle
import gzip
#from sparray import sparray


### PARSE ONE STAGE ONLY!!! (too much memory)
# main_dir = '/Users/kasia/Documents/PhD/outputs/Zv9/wiggle/0_2-4Cell'
main_dir = '/Users/kasia/Documents/PhD/outputs/Zv9/WIG_split_by_length/5_Bud'
# 0_2-4Cell     X
# 1_256Cell     X
# 2_1KCell      X
# 3_Dome        split!!
# 4_Shield      split!!
# 5_Bud         X
# 6_28hpf       split!!
# 7_5dpf        split!!


# just paths to files here
fwd_handles = []
for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(main_dir, 'fwd')))):
	fwd_handles.extend([ParseWig(os.path.join(main_dir, 'fwd', file))])

rev_handles = []
for file in (filter(lambda f: not f.startswith('.'), os.listdir(os.path.join(main_dir, 'rev')))):
    rev_handles.extend([ParseWig(os.path.join(main_dir, 'rev', file))])




transcripts_lengths = {}

for tr in orfs:
    strand = geneToInterval[tr][0].strand
    if strand == '+':
        handles = fwd_handles
    elif strand == '-':
        handles = rev_handles
    ### initiate sequence array
    sequencearray = np.zeros((len(fwd_handles), geneLength[tr])) # IF USING FWD_HANDLES!
    for i in range(len(handles)): # for different read lengths
        cds = np.array([])
        for interval in geneToInterval[tr]:
            exonint = np.asarray(handles[i].fetch_all_scores(interval.chrom, interval.start+1, interval.end+1))
            cds = np.append(cds, exonint)
            #put into sequencearray
        cds[np.isnan(cds)] = 0 # get all exons for one handle
        sequencearray[i] = cds ###lengths???
    ### get ORF positions
    orfarray = sequencearray[:, orfs[tr][0][0]:orfs[tr][0][1]]
    # for + strand, sparse
    sp_orf = sparse.csr_matrix(orfarray)
    sp_tran = sparse.csr_matrix(sequencearray)
    if strand == '-':
        orfarray = np.fliplr(orfarray)
        sequencearray = np.fliplr(sequencearray)
        sp_orf = sparse.csr_matrix(orfarray)
        sp_tran = sparse.csr_matrix(sequencearray)
    ### put all into a dictionary
    transcripts_lengths[tr] = sp_tran, sp_orf



pickle.dump(transcripts_lengths, gzip.open('../meta_data/split_by_length/sparse/5_Bud.p.gz', 'wb'))
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################


orfdict = {k: np.zeros((8, 7, orfs[k][0][1]-orfs[k][0][0])) for k in tbl.keys()}

for key in tbl:
    for i in range(len(tbl[key])):
        for j in range(len(tbl[key][i])):
            orf = tbl[key][i][j, orfs[key][0][0]:orfs[key][0][1]]
            if (geneToInterval[key][0].strand == '-'):
                orf = orf[::-1]
            orfdict[key][i][j] = orf


# merge stages
# instead of having >> 8_stages x 7_lengths x geneLength << for each key
# have 7 x geneLength
merged_by_length = {}
for key in orfdict:
    lengths = np.zeros((7, orfs[key][0][1]-orfs[key][0][0]))
    for i in range(len(orfdict[key])):
        for j in range(len(orfdict[key][i])):
            lengths[j] += orfdict[key][i][j]
    merged_by_length[key] = lengths


# calculate codon coverage
codon_cov = {}
for key in merged_by_length:
    m = merged_by_length[key]
    if (len(m[0]) % 3) == 0:
        pep_array = np.zeros((7, len(m[0])/3))
        for i in range(len(m)):
            for j in range(0, len(m[i]), 3):
                scc = abs(m[i][j] + m[i][j+1] + m[i][j+2])
                pep_array[i][j/3] = scc
            codon_cov[key] = pep_array

# merge the lengths, to calculate frequencies and choose stall sites
codon_cov_merged = {}
for key in codon_cov:
    merged = np.zeros(len(codon_cov[key][0]))
    for i in range(len(codon_cov[key])):
        merged += codon_cov[key][i]
    codon_cov_merged[key] = merged

# relative frequencies at given codon
freqs = {}
for key in codon_cov_merged:
    freq_array = np.zeros(len(codon_cov[key][0]))
    length = len(codon_cov_merged[key])
    sumcov = sum(codon_cov_merged[key])
    for i in range(len(codon_cov_merged[key])):
        cov = codon_cov_merged[key][i]
        #f = (cov * (1.0/length)) / (float(sumcov)/length) #it's just cov / sumcov ...
        f = (float(cov)/float(sumcov)) / length
        freq_array[i] = f
    freqs[key] = freq_array


# choose frequency cutoff
#threshold = 0.05
threshold = 0.0002

# save stall positions - above the threshold
genes_with_peaks = {}
for key in freqs:
    if max(freqs[key]) > threshold:
        peak_positions = []
        for i in range(len(freqs[key])):
            if freqs[key][i] > threshold:
                peak_positions += [i*3] ### saves as 1st nucleotide in codon, distance from ORF start
        genes_with_peaks[key] = peak_positions


# extract -30 +33 windows over stall positions
# if so long are not available - fill with zeros
# do it for 7 separate lengths (25-31) - for later plotting
lengths = {}
for i in range (25, 32):
    lengths[i] = np.array([0]*63)

for key in genes_with_peaks:
    for peak in genes_with_peaks[key]:
        start = peak-30
        end = peak+33
        for l in range(7):
            if start < 0:
                #lengths[l+25][0:abs(start)] = 0
                lengths[l+25][abs(start):] += merged_by_length[key][l][0:end]
            elif end > len(merged_by_length[key][l]):
                lengths[l+25][0:63-(end-len(merged_by_length[key][l]))] += merged_by_length[key][l][start:len(merged_by_length[key][l])]
                #lengths[l+25][len(merged_by_length[key][l]):-1] = 0
            else:
                lengths[l+25] += merged_by_length[key][l][start:end]



# plot stacked histogram
colors = ['#660000', '#CC0000', '#FF6666', '#FFFFFF', '#CCCCCC', '#666666', 'k']
x = range(-30,33)
p1 = plt.bar(x, stall[25], color=colors[0])
p2 = plt.bar(x, stall[26], color=colors[1])
p3 = plt.bar(x, stall[27], color=colors[2])
p4 = plt.bar(x, stall[28], color=colors[3])
p5 = plt.bar(x, stall[29], color=colors[4])
p6 = plt.bar(x, stall[30], color=colors[5])
p7 = plt.bar(x, stall[31], color=colors[6])
plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0]), ('25', '26', '27', '28', '29', '30', '31'))



# sum total coverage per certain length
mbl = pickle.load(gzip.open('merged_by_length.p.gz', 'rb'))

summed = np.array([0]*7)
for key in mbl:
    for i in range(len(mbl[key])):
        summed[i] += sum(mbl[key][i])

### RESULT: array([ 1481129,  4678173, 10118275, 13143870,  9584410,  4401804,  1380935])
plt.bar(range(25, 32), summed, align='center')