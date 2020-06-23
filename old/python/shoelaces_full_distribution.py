#
# INPUT: A FASTA and BED file of a gene set and the reads as BAM
#
# OUTPUT: 
#
#
#
import argparse, sys, os
import HTSeq
from HTSeq import GenomicPosition, GenomicInterval
import pysam
from subprocess import Popen, PIPE

from math import log
import numpy



def LENGTHscore(real, noise):
#    sys.stderr.write("R(%s) N(%s)\n" % (sum(real), sum(noise)))

    realReads  = sum(real)
    noiseReads = sum(noise)
    return log(float(realReads+1) / float(noiseReads+1), 2)


def ORFscore(mergedProfile, tot):
    frame1 = [mergedProfile[i] for i in range(0,len(mergedProfile),3)]
    frame2 = [mergedProfile[i] for i in range(1,len(mergedProfile),3)]
    frame3 = [mergedProfile[i] for i in range(2,len(mergedProfile),3)]

    tot3 = tot/3
    
    sf1 = sum(frame1)
    f1 = (sf1 - tot3)
    f1 *= f1

    sf2 = sum(frame2)
    f2 = (sf2 - tot3)
    f2 *= f2

    sf3 = sum(frame3)
    f3 = (sf3 - tot3)
    f3 *= f3

    if tot3 == 0:
        score = 0
    else:
        score = log( (f1+f2+f3)/tot3 + 1, 2)
        ## FIX:CHECK +1 SCORE

    if sf1 < sf2 or sf1 < sf3:
        score = -score

    return (score, sum(frame1), sum(frame2), sum(frame3))



def summarize_ORFs(name, orfs, profile, strand, allowed_lengths, transcript_length, exp, out, distOutput):

    for orf in orfs:
        dist = {}
        last_dist = {}
        outdist = {}

        sumlen = {}
        cvg = numpy.zeros( orf[1]-orf[0], dtype='i' )
        for bamfile in profile:
            for length in profile[bamfile]:       
                sumlen[length] = sum(profile[bamfile][length][orf[0]:orf[1]])
                if length in allowed_lengths:
                    cvg += profile[bamfile][length][orf[0]:orf[1]]

        orfcvg = cvg[orf[0]:orf[1]]
        upstream = cvg[:orf[0]]
        downstream = cvg[orf[1]:]

        upstream_tot = sum(upstream)
        downstream_tot = sum(downstream)
        
        inside_tot = sum(orfcvg)

        orf_score, f1, f2, f3 = ORFscore(orfcvg, inside_tot)
        
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, orf[0], orf[1], transcript_length, inside_tot, upstream_tot, downstream_tot, f1, f2, f3, orf_score))

        ky = sorted(sumlen.keys())
        u = ["%s" % (str(sumlen[x])) for x in ky]        
        distOutput.write("%s\t%s\t%s\t%s\n" % (name, orf[0], orf[1], "\t".join(u)))

def summarize_genes(gene_chains, chainLength, orfs, lengths, plus_wig, minus_wig, outputDir):
    allowed_lengths = ['27','28','29','30']
    exp = None

    out = open(outputDir + "/orfs.table", "w")
    dist = open(outputDir + "/length_distribution.table", "w")

    for name, chain in gene_chains.iteritems():
        offset = 0
        added = 0

        # Skip transcripts without CDS
        if len(chain) == 0:
            continue

        strand = chain[0].strand

        coverage = {}
        coverage["TEST"] = {}

        for length in lengths:
            coverage["TEST"][length] = numpy.zeros( chainLength[name], dtype='i' )

        for length in lengths:
            current = 0
            for interval in chain:
                if strand == "+":
                    wigfile = plus_wig[length]
                elif strand == "-":
                    wigfile = minus_wig[length]


#                sys.stderr.write("bigWigSummary %s %s %s %s %s \n" % (wigfile, interval.chrom, interval.start, interval.end, interval.end-interval.start))
                prog = Popen("bigWigSummary %s %s %s %s %s " % (wigfile, interval.chrom, interval.start, interval.end, interval.end-interval.start), shell=True, stdout=PIPE)    
                output = prog.communicate()  

                if (prog.returncode != 0):
#                sys.stderr.write("Interval failed: %s\n" % interval);
                    continue
            
                output = output[0]
                output = output.rstrip()
                o = map (lambda x : 0 if x == "n/a" else float(x), output.split("\t"))
            
#                sys.stderr.write("\tCL(%s) %s - %s  O(%s)\n" % (len(coverage["TEST"][length]), current, current+len(o), len(o)))
#            sys.stderr.write("O: %s\n" % o)
                coverage["TEST"][length][current:(current+len(o))] = o
                current += len(o)

#        sys.stderr.write("%s: %s\n" % (name, len(coverage)))
#       sys.stderr.write("COV: %s\n\n" % coverage)

        if orfs.has_key(name):
            summarize_ORFs(name, orfs[name], coverage, gene_chains[name][0].strand, allowed_lengths, chainLength[name], exp, out, dist)
        else:
            sys.stderr.write("No ORFs: %s\n" % name)

    out.close()
    dist.close()

def make_shoelace_stats(bed_file, lengths, plus_wig, minus_wig, outputDir):
    intervals, geneToInterval, geneLength = BED_to_intervals(args.BED)    
    orfs = BED_get_CDS_mRNA_coord(args.BED)   
    summarize_genes(geneToInterval, geneLength, orfs, lengths, plus_wig, minus_wig, outputDir)



def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("BED", help="BED coordinates of the transcripts", metavar="BED")
    parser.add_argument("LENGTHS", help="Lengths of reads", metavar="LENGTHS")
    parser.add_argument("BW_PLUS", help="WIG file with reads starts on plus strand", metavar="BW_PLUS")
    parser.add_argument("BW_MINUS", help="WIG file with reads starts on minus strand", metavar="BW_MINUS")
    parser.add_argument("-o", "--outputDir", default=".",  help="Output directory")

    return parser.parse_args()


def BED_to_intervals(bed_file):
    bed = open(bed_file, "rU")

    intervals = set()
    geneToInterval = {}
    geneLength = {}

    for line in bed:
        gene = line.split("\t")
        chrom = gene[0]
        start = int(gene[1])
        strand = gene[5]
        
        geneStr = gene[3] #"%s:::%s" % (gene[3], gene[1])
        if geneToInterval.has_key(geneStr):
            sys.stderr.write("Warning: duplicate gene names: %s. Overwriting previous entry.\n" % geneStr)
        geneToInterval[geneStr] = []
        geneLength[geneStr] = 0

        exonLengths = gene[10].split(",")
        exonStarts = gene[11].split(",")

        if gene[10][-1] == ",":
            exonLengths.pop()
        if gene[11][-1] == ",":
            exonStarts.pop()

        offset = 0
        for exonStart, exonLen in zip(exonStarts, exonLengths):
            s = start + int(exonStart)
            e = s + int(exonLen)

            interval = GenomicInterval(chrom, s, e, strand)
            intervals.add(interval)

            geneToInterval[geneStr].append(interval)
            geneLength[geneStr] += interval.length

    return intervals, geneToInterval, geneLength

def BED_get_CDS_mRNA_coord(bed_file):
    bed = open(bed_file, "rU")
    cds = {}

    # Make unique set of CDSs
    for line in bed:
        gene = line.split("\t")
        strand = gene[5]
#        geneStr = "%s:::%s" % (gene[3], gene[1])
        geneStr = gene[3]

        txStart = int(gene[1])
        cds_start = int(gene[6]) - txStart
        cds_end = int(gene[7]) - txStart

        # Non-coding
        if cds_start == cds_end:
            continue

        exonLengths = gene[10].split(",")
        exonStarts = gene[11].split(",")

        # BED files occasionally end with a ","
        if gene[10][-1] == ",":
            exonLengths.pop()
        if gene[11][-1] == ",":
            exonStarts.pop()

        offset = 0
        for exonStart, exonLen in zip(exonStarts, exonLengths):
            exonStart = int(exonStart)
            exonLen = int(exonLen)

            if cds_start >= exonStart and cds_start <= exonStart+exonLen:
                new_cds_start = offset + (cds_start - exonStart)

            if cds_end >= exonStart and cds_end <= exonStart+exonLen:
                new_cds_end = offset + (cds_end - exonStart)

            offset += exonLen

        if new_cds_start > new_cds_end:
            sys.stderr.write("CDS length is negative for transcript %s\n" % gene[3])
            continue

        if not cds.has_key(geneStr):
            cds[geneStr] = []
        cds[geneStr].append((new_cds_start, new_cds_end))

    return cds


def BED_get_CDS_genome_coord(bed_file):
#    sys.stderr.write("Parsing gene annotation: %s\n" % bed_file)

    bed = open(bed_file, "rU")
    startpos = set()
    endpos = set()

    # Make unique set of CDSs
    for line in bed:
        gene = line.split("\t")
        strand = gene[5]

        if strand == "+":
            start = gene[6]
            end = gene[7]
        else :
            start = gene[7]
            end = gene[6]
        
        startpos.add(GenomicPosition(gene[0], int(start), strand))
        endpos.add(GenomicPosition(gene[0], int(end), strand))

    return startpos, endpos

def main(args):
    exp_file = None

    plus_bw = {}
    minus_bw = {}
    lengths = args.LENGTHS.split(",")
    for f,l in zip(args.BW_PLUS.split(","), lengths):
        plus_bw[l] = f
    for f,l in zip(args.BW_MINUS.split(","), lengths):
        minus_bw[l] = f


    make_shoelace_stats(args.BED, lengths, plus_bw, minus_bw, args.outputDir)    


if __name__ == '__main__':
     args = parse_options()

     if not os.path.exists(args.outputDir):
         os.makedirs(args.outputDir)

     main(args)

