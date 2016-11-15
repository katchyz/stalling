from HTSeq import GenomicPosition, GenomicInterval
import cPickle as pickle
import gzip
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord, BED_to_intervals
from scipy import sparse
import numpy as np
import os
from wiggle import ParseWig


bedpath = '/Home/ii/katchyz/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
gtfpath = '/Home/ii/katchyz/DATA/genomes/GTF/Mus_musculus.GRCm38.79.chr.gtf'

path_to_feature_coord = '/Home/ii/katchyz/META/mouse_leader_CDS_trailer_coord.p.gz'


stage_dirs = ['EB_GA_1', 'EB_GA_2', 'EB_GA_3', 'EB_HiSeq', 'ESC_36h_GA_1', 'ESC_36h_GA_2', 'ESC_36h_GA_3', \
'ESC_36h_HiSeq', 'ESC_ff_GA_1', 'ESC_ff_GA_2', 'ESC_ff_GA_3', 'ESC_ff_HiSeq', 'ESC_ff_chx_GA_1', 'ESC_ff_chx_GA_2', \
'ESC_ff_emet_GA', 'mouse_none']


wigs_dir = '/Home/ii/katchyz/OUT/mouse_out/by_length/'

fr_lengths = range(25, 36)

lct_dict_path_dump = '/Home/ii/katchyz/META/leader_cds_trailer/mouse.p.gz'

counts_path_dump = '/Home/ii/katchyz/META/leader_cds_trailer/counts_mouse.p.gz'


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


startpos, endpos = BED_get_CDS_genome_coord(bedpath)

orfs = BED_get_CDS_mRNA_coord(bedpath)
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)


# read gtf, get a set with gene names and corresponding transcript names
otpg = {}
gtf = open(gtfpath)

for line in gtf:
    if not line.startswith('#') and line.split()[2] != 'gene':
        gene = line.split()[9][1:-2]
        tr = line.split()[13][1:-2]
        # print gene, tr
        if gene not in otpg:
            otpg[gene] = [tr]
        elif gene in otpg and tr not in otpg[gene]:
            otpg[gene].append(tr)


# for each of the genes, merge transcript intervals to get assigned CDSs, remaining regions allocate to leaders and trailers
# (save coordinates)

def merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged



feature_coord = {}

for gene in otpg:
    l, c, t = [], [], []
    for tr in otpg[gene]:
        if geneToInterval[tr][0].strand == '+':
            startF, stopF = 0, 0
            for i in geneToInterval[tr]:
                if i.end > startpos[tr].start and startF == 0: # exon with start codon
                    startF = 1
                    c.append((startpos[tr].start, i.end))
                    l.append((i.start, startpos[tr].start))
                elif i.end >= endpos[tr].end and startF == 1: # exon with stop codon
                    stopF = 1
                    c.append((i.start, endpos[tr].end))
                    t.append((endpos[tr].end, i.end))
                else:
                    if startF == 1 and stopF == 0: # fully coding exons
                        c.append((i.start, i.end))
                    elif startF == 0 and stopF == 0: # leader
                        l.append((i.start, i.end))
                    elif startF == 1 and stopF == 1: # trailer
                        t.append((i.start, i.end))
        elif geneToInterval[tr][0].strand == '-':
            startF, stopF = 0, 0
            for i in geneToInterval[tr]:
                if i.end > endpos[tr].start and stopF == 0: # exon with stop codon
                    stopF = 1
                    c.append((endpos[tr].start, i.end))
                    t.append((i.start, endpos[tr].start))
                elif i.end >= startpos[tr].end and stopF == 1: # exon with start codon
                    startF = 1
                    c.append((i.start, startpos[tr].end))
                    l.append((startpos[tr].end, i.end))
                else:
                    if startF == 0 and stopF == 1: # fully coding exons
                        c.append((i.start, i.end))
                    elif startF == 0 and stopF == 0: # trailer
                        t.append((i.start, i.end))
                    elif startF == 1 and stopF == 0: #leader
                        l.append((i.start, i.end))
    #print 'leader: ', merge_intervals(l), 'CDS: ', merge_intervals(c), 'trailer: ', merge_intervals(t)
    feature_coord[gene] = merge_intervals(l), merge_intervals(c), merge_intervals(t), geneToInterval[tr][0].strand, geneToInterval[tr][0].chrom



features = {}
for gene in feature_coord:
    l = feature_coord[gene][0]
    c = feature_coord[gene][1]
    t = feature_coord[gene][2]
    strand = feature_coord[gene][3]
    chrom = feature_coord[gene][4]
    leader, trailer = [], []
    if strand == '+':
        CDSstart = c[0][0]
        CDSend = c[-1][1]
        for i in l:
            if CDSstart >= i[0] and CDSstart <= i[1]:
                leader.append((i[0], CDSstart))
                break
            else:
                leader.append((i[0], i[1]))
        for j in t:
            if CDSend >= j[0] and CDSend <= j[1]:
                trailer.append((CDSend, j[1]))
            elif CDSend < j[0]:
                trailer.append((j[0], j[1]))
    elif strand == '-':
        CDSstart = c[-1][1]
        CDSend = c[0][0]
        for i in t:
            if CDSend >= i[0] and CDSend <= i[1]:
                trailer.append((i[0], CDSend))
                break
            else:
                trailer.append((i[0], i[1]))
        for j in l:
            if CDSstart >= j[0] and CDSstart <= i[1]:
                leader.append((CDSstart, j[1]))
            elif CDSstart < j[0]:
                leader.append((j[0], j[1]))
    print gene, 'leader: ', leader, 'CDS: ', c, 'trailer: ', trailer, 'strand: ', strand, 'chrom: ', chrom
    features[gene] = leader, c, trailer, strand, chrom


pickle.dump(features, gzip.open(path_to_feature_coord, 'wb'))

#################################
#################################
#################################

# on actual data:
# extract leader, cds and trailer regions from WIGs

def extract_intervals_from_wigs(dirpath, features, filename, fragment_length):
    """ Dumps transcript coverage into python dictionary """
    fl = str(fragment_length)
    fwd_handle = ParseWig(os.path.join(dirpath, 'fwd', filename+'-'+fl+'-forward.wig'))
    rev_handle = ParseWig(os.path.join(dirpath, 'rev', filename+'-'+fl+'-reverse.wig'))
    gene_features = {}
    for gene in features:
        strand = features[gene][3]
        if strand == '+':
            handle = fwd_handle
        elif strand == '-':
            handle = rev_handle
        ### initiate sequence array
        leader = np.array([])
        cds = np.array([])
        trailer = np.array([])
        if features[gene][0]:
            for i in features[gene][0]: #leaders
                exonint = np.asarray(handle.fetch_all_scores(features[gene][4], i[0]+1, i[1]+1))
                leader = np.append(leader, exonint)
            leader[np.isnan(leader)] = 0
        if features[gene][1]:
            for i in features[gene][1]: #CDSs
                exonint = np.asarray(handle.fetch_all_scores(features[gene][4], i[0]+1, i[1]+1))
                cds = np.append(cds, exonint)
            cds[np.isnan(cds)] = 0
        if features[gene][2]:
            for i in features[gene][2]: #trailers
                exonint = np.asarray(handle.fetch_all_scores(features[gene][4], i[0]+1, i[1]+1))
                trailer = np.append(trailer, exonint)
            trailer[np.isnan(trailer)] = 0
        # minus strand
        if strand == '-':
            leader = leader[::-1]
            cds = cds[::-1]
            trailer = trailer[::-1]
        if leader.any():
            sp_leader = sparse.csr_matrix(leader)
        else:
            sp_leader = sparse.csr_matrix([0])
        if cds.any():        
            sp_cds = sparse.csr_matrix(cds)
        else:
            sp_cds = sparse.csr_matrix([0])
        if trailer.any():
            sp_trailer = sparse.csr_matrix(trailer)
        else:
            sp_trailer = sparse.csr_matrix([0])
        gene_features[gene] = sp_leader, sp_cds, sp_trailer
    return gene_features



# ...OR: get it from transcript dictionary? (not that easy)
# save a dictionary for each gene, with leader, cds, trailer (for different fragment length)

d = {}
for stage in stage_dirs:
    d[stage] = {}
    for fr_len in fr_lengths:
        main_dir = wigs_dir+stage
        d[stage][str(fr_len)] = extract_intervals_from_wigs(main_dir, features, stage, fr_len)
        print 'read ', stage, fr_len, '\n'


### merge different lengths

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
    for tr in d[stage][str(fr_lengths[0])]:
        seq_orf[stage][tr] = d[stage][str(fr_lengths[0])][tr][0], d[stage][str(fr_lengths[0])][tr][1], d[stage][str(fr_lengths[0])][tr][2]

for stage in seq_orf:
    for l in fr_lengths[1::]:
        for tr in d[stage][str(l)]:
            csr_vappend(seq_orf[stage][tr][0], d[stage][str(l)][tr][0])
            csr_vappend(seq_orf[stage][tr][1], d[stage][str(l)][tr][1])
            csr_vappend(seq_orf[stage][tr][2], d[stage][str(l)][tr][2])


pickle.dump(seq_orf, gzip.open(lct_dict_path_dump, 'wb'))


# no_fl = len(fr_lengths)

# counts = {}
# for stage in seq_orf:
#     counts[stage] = np.zeros(no_fl), np.zeros(no_fl), np.zeros(no_fl)
#     for tr in seq_orf[stage]:
#         if seq_orf[stage][tr][0].shape[1] > 0:
#             for l in range(no_fl):
#                 counts[stage][0][l] += seq_orf[stage][tr][0][l].sum() # leader
#         for l in range(no_fl):
#             counts[stage][1][l] += seq_orf[stage][tr][1][l].sum() # CDS
#         if seq_orf[stage][tr][2].shape[1] > 0:
#             for l in range(no_fl):
#                 counts[stage][2][l] += seq_orf[stage][tr][2][l].sum() # trailer


# pickle.dump(counts, gzip.open(counts_path_dump, 'wb'))

###
# save data frames to make plots
# from pandas import DataFrame

# for stage in counts:
#     data = np.vstack((counts[stage][0], counts[stage][1], counts[stage][2]))
#     df = DataFrame(data=data)
#     df.to_csv('/Volumes/USELESS/META/leader_cds_trailer/'+stage+'.csv')




