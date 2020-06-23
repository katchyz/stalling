### STALLING - human
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

### import (separated by length?)
# andreev2015 "Andreev.+"
# fritsch2012 "Fritsch.+"
# gonzalez2014 "Gonzalez.+"
# guo2010 "Guo.+"
# hsieh2012 "Hsieh.+"
# ingolia2012 "Ingolia_NT_2012.+"
# ingolia2014 "Ingolia_NT_2014.+"
# lee2012
# liu2013
# rutkowski2015
# sidrauski2015
# stern2012
# stumpf2013
# subtelny2014
# yoon2014

### GTF
txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# strand
strand <- transcripts(txdb)
strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
# exons
exons <- exonsBy(txdb, by="tx", use.names=TRUE)


####################################################
# make empty granges
starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                                  IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                                  Rle(txLengths$strand, txLengths$tx_len))

empty_gr$riboseq <- rep(0, length(empty_gr))


########################################################################################################

#pattern <- "Andreev.+"
#save_prefix <- "andreev2015_"

pattern <- "Gonzalez.+"
save_prefix <- "gonzalez2014_"

wig_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/tracks_per_length"
libs <- list.files(path = wig_path, pattern = pattern)

srr <- sapply(libs, function(x){substr(x,(nchar(x)-24),(nchar(x)-15))})
lengths <- sapply(libs, function(x){substr(x,(nchar(x)-13),(nchar(x)-12))})
srr_len <- data.frame(srr = srr, lengths = lengths)

runs <- unique(srr)
for (i in 1:length(runs)) {
  run = runs[i]
  lens = as.character(unique(srr_len[srr_len$srr == run,]$lengths))
  ## get empty granges
  gr <- empty_gr
  for (j in 1:length(lens)) {
    l = lens[j]
    # load library
    # fw
    pattern_fw <- paste0(pattern, run, "-", l, "-forward.wig")
    wig_fw <- list.files(path = wig_path, pattern = pattern_fw)
    wig_fw <- import.wig(c(file.path(wig_path, wig_fw)))
    strand(wig_fw) <- rep("+", length(wig_fw))
    # rv
    pattern_rv <- paste0(pattern, run, "-", l, "-reverse.wig")
    wig_rv <- list.files(path = wig_path, pattern = pattern_rv)
    wig_rv <- import.wig(c(file.path(wig_path, wig_rv)))
    strand(wig_rv) <- rep("+", length(wig_rv))
    # join
    wig <- c(wig_fw, wig_rv)
    # map to tx
    wig_tx <- mapToTranscripts(wig, exons, ignore.strand = FALSE)
    mcols(wig_tx)[[l]] <- wig[wig_tx$xHits]$score
    # find overlaps
    fo <- findOverlaps(gr, wig_tx, type="equal")
    lenL <- rep(0, length(gr))
    lenL[from(fo)] <- mcols(wig_tx)[[l]][to(fo)]
    mcols(gr)[[l]] <- lenL
    # add length to total
    gr$riboseq <- gr$riboseq + lenL
  }
  ## save with run name
  save_file <- paste0("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges/", save_prefix, run, ".Rsave")
  save(gr, file = save_file)
}



#wig_tx <- split(wig_tx, seqnames(wig_tx))
#wig_tx <- wig_tx[order(names(wig_tx))]


