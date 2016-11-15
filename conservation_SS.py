 ### GET STALL SITE CONSERVATION
import cPickle as pickle
import gzip
from bx.bbi.bigwig_file import BigWigFile
from shoelaces_full_distribution import BED_to_intervals, BED_get_CDS_mRNA_coord
from HTSeq import GenomicPosition, GenomicInterval
import itertools
import numpy as np

# get mouse (conserved) stall site positions (and all too...)
css = pickle.load(gzip.open('/Volumes/USELESS/META/conserved_stall_sites.p.gz', 'rb'))
gwp = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_mouse.p.gz', 'rb'))
# get their genomic coordinates, chromosome # (strand?)
bed = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'

phastcons = '/Volumes/USELESS/DATA/mm10.60way.phastCons.bw'
bwh = BigWigFile(open(phastcons, "rb"))

phylop = '/Volumes/USELESS/DATA/mm10.60way.phyloP60way.bw'
bwh = BigWigFile(open(phylop, "rb"))

def BED_get_CDS_genome_coord(bed_file):
#    sys.stderr.write("Parsing gene annotation: %s\n" % bed_file)
    bed = open(bed_file, "rU")
    startpos = {}
    endpos = {}
    # Make unique set of CDSs
    for line in bed:
        gene = line.split("\t")
        strand = gene[5]
        geneStr = gene[3]
        if strand == "+":
            start = gene[6]
            end = gene[7]
        else :
            start = gene[7]
            end = gene[6]
        startpos[geneStr] = GenomicPosition(gene[0], int(start), strand)
        endpos[geneStr] = GenomicPosition(gene[0], int(end), strand)
    return startpos, endpos

intervals, geneToInterval, geneLength = BED_to_intervals(bed)
mRNA_coord = BED_get_CDS_mRNA_coord(bed)
startpos, endpos = BED_get_CDS_genome_coord(bed)


SS_CDS_coord = {}
for gene in css:
	for peak in css[gene]:
		if 'mouse' in css[gene][peak]:
			# get transcript coordinates
			# print gene, css[gene][peak]['mouse']
			tr = css[gene][peak]['mouse'][1]
			if not tr in SS_CDS_coord:
				SS_CDS_coord[tr] = [css[gene][peak]['mouse'][0]]
			else:
				SS_CDS_coord[tr].append(css[gene][peak]['mouse'][0])


# converts list of integers to ranges
def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]


SS_cons = {}
for tr in SS_CDS_coord:
	SS_cons[tr] = {}
	strand = geneToInterval[tr][0].strand
	chrom = geneToInterval[tr][0].chrom
	tr_gen_int = []
	for interval in geneToInterval[tr]:
		tr_gen_int.extend(range(interval.start, interval.end))
	# get the right position for + and minus strand
	if strand == '+':
		for peak in SS_CDS_coord[tr]:
			SS_cons[tr][peak] = np.array([])
			pos_on_tr = int(peak) + mRNA_coord[tr][0][0]
			# get this +/-15
			rans = ranges(tr_gen_int[pos_on_tr-15:pos_on_tr+15])
			for ran in rans:
				SS_cons[tr][peak] = np.append(SS_cons[tr][peak], bwh.get_as_array(chrom, ran[0], ran[-1]))
	if strand == '-':
		for peak in SS_CDS_coord[tr]:
			SS_cons[tr][peak] = np.array([])
			pos_on_tr = -int(peak) + geneLength[tr] - mRNA_coord[tr][0][1] -3 # starting from the end
			# get this +/-15
			rans = ranges(tr_gen_int[-pos_on_tr-15:-pos_on_tr+15])
			for ran in rans:
				SS_cons[tr][peak] = np.append(SS_cons[tr][peak], bwh.get_as_array(chrom, ran[0], ran[-1]))


for tr in SS_cons:
	for peak in SS_cons[tr]:
		if SS_cons[tr][peak].any():
			SS_cons[tr][peak][np.isnan(SS_cons[tr][peak])] = 0


####### CONTROL
# get control - positions 50nt downstream
SS_cont = {}
for tr in SS_CDS_coord:
	SS_cont[tr] = {}
	strand = geneToInterval[tr][0].strand
	chrom = geneToInterval[tr][0].chrom
	tr_gen_int = []
	for interval in geneToInterval[tr]:
		tr_gen_int.extend(range(interval.start, interval.end))
	# get the right position for + and minus strand
	if strand == '+':
		for peak in SS_CDS_coord[tr]:
			SS_cont[tr][peak] = np.array([])
			pos_on_tr = int(peak) + mRNA_coord[tr][0][0] + 50
			# get this +/-15
			rans = ranges(tr_gen_int[pos_on_tr-15:pos_on_tr+15])
			for ran in rans:
				SS_cont[tr][peak] = np.append(SS_cont[tr][peak], bwh.get_as_array(chrom, ran[0], ran[-1]))
	if strand == '-':
		for peak in SS_CDS_coord[tr]:
			SS_cont[tr][peak] = np.array([])
			pos_on_tr = -int(peak) + geneLength[tr] - mRNA_coord[tr][0][1] -3 - 50# starting from the end
			# get this +/-15
			rans = ranges(tr_gen_int[-pos_on_tr-15:-pos_on_tr+15])
			for ran in rans:
				SS_cont[tr][peak] = np.append(SS_cont[tr][peak], bwh.get_as_array(chrom, ran[0], ran[-1]))


for tr in SS_cont:
	for peak in SS_cont[tr]:
		if SS_cont[tr][peak].all():
			SS_cont[tr][peak][np.isnan(SS_cont[tr][peak])] = 0

#################

means = []
for tr in SS_cons:
	for peak in SS_cons[tr]:
		if SS_cons[tr][peak].any():
			means.append(np.mean(SS_cons[tr][peak]))

ends = np.array([])
for tr in SS_cons:
	strand = geneToInterval[tr][0].strand
	for peak in SS_cons[tr]:
		if SS_cons[tr][peak].all():
			if strand == '+':
				ends = np.append(ends, SS_cons[tr][peak][0:10])
			elif strand == '-':
				ends = np.append(ends, SS_cons[tr][peak][-10:-1])

cends = np.array([])
for tr in SS_cont:
	strand = geneToInterval[tr][0].strand
	for peak in SS_cont[tr]:
		if SS_cont[tr][peak].all():
			if strand == '+':
				cends = np.append(cends, SS_cont[tr][peak][0:10])
			elif strand == '-':
				cends = np.append(cends, SS_cont[tr][peak][-10:-1])



for tr in SS_cont:
	strand = geneToInterval[tr][0].strand
	for peak in SS_cont[tr]:
		if SS_cont[tr][peak].all():
			if strand == '+':
				print SS_cont[tr][peak].tolist()


import matplotlib.pyplot as plt
plt.hist(means)

cmeans = []
for tr in SS_cont:
	print tr
	for peak in SS_cont[tr]:
		print peak
		if SS_cont[tr][peak].all():
			cmeans.append(np.mean(SS_cont[tr][peak]))

for i in range(len(cmeans)):
	if np.isnan(cmeans[i]):
		cmeans[i] = 0

plt.hist(cmeans)


# bwh.get_as_array('chr2', 180971999, 180972028)
# bwh.get('chr1', 1000000, 1000100)

# extract from BigWig file
# (BW is 1-based, my python SSs are 0-based)
# extract regions to compare it with (exonic, 50nt upstream or sth)

# compare...



