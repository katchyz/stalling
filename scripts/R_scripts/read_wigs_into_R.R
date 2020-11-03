## read in riboseq wiggle files, put different lengths into granges

library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"

##### YEAST #####
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Saccharomyces_cerevisiae.R64-1-1.79.gtf")), format="gtf")
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


# make empty granges
starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                    IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                    Rle(txLengths$strand, txLengths$tx_len))

empty_gr$riboseq <- rep(0, length(empty_gr))



save_gr <- function(empty_gr, pattern, wig_path, exons, save_prefix) {
  libs <- list.files(path = wig_path, pattern = pattern)
  lengths <- sapply(libs, function(x){substr(x,(nchar(x)-13),(nchar(x)-12))})
  lens = unique(lengths)
  ## get empty granges
  gr <- empty_gr
  for (j in 1:length(lens)) {
    l = lens[j]
    # load library
    # fw
    pattern_fw <- paste0(pattern, "-", l, "-forward.wig")
    wig_fw <- list.files(path = wig_path, pattern = pattern_fw)
    wig_fw <- import.wig(c(file.path(wig_path, wig_fw)))
    strand(wig_fw) <- rep("+", length(wig_fw))
    # rv
    pattern_rv <- paste0(pattern, "-", l, "-reverse.wig")
    wig_rv <- list.files(path = wig_path, pattern = pattern_rv)
    wig_rv <- import.wig(c(file.path(wig_path, wig_rv)))
    strand(wig_rv) <- rep("-", length(wig_rv))
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
  save_file <- paste0("/Volumes/USELESS/STALLING/gr/", save_prefix, ".Rsave")
  save(gr, file = save_file)
}

##########
wig_path <- "/Volumes/USELESS/STALLING/wigs/yeast"
pattern <- "gb15_meio1"
save_prefix <- "yeast_gb15_meio1"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "gb15_meio2"
save_prefix <- "yeast_gb15_meio2"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "yeast_chx"
save_prefix <- "yeast_chx"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "yeast_none"
save_prefix <- "yeast_none"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "yeast_subtelny"
save_prefix <- "yeast_subtelny"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)
################################################################################################################

##### FRUIT FLY #####
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Drosophila_melanogaster.BDGP6.79.gtf")), format="gtf")
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

# make empty granges
starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                    IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                    Rle(txLengths$strand, txLengths$tx_len))

empty_gr$riboseq <- rep(0, length(empty_gr))

##########
wig_path <- "/Volumes/USELESS/STALLING/wigs/fruitfly"
pattern <- "0-2h_embryo_A"
save_prefix <- "fruitfly_0-2h_embryo_A"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "0-2h_embryo_B"
save_prefix <- "fruitfly_0-2h_embryo_B"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "S2cell_150_A"
save_prefix <- "fruitfly_S2cell_150_A"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "S2cell_150_B"
save_prefix <- "fruitfly_S2cell_150_B"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "S2cell_250_A"
save_prefix <- "fruitfly_S2cell_250_A"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "S2cell_250_B"
save_prefix <- "fruitfly_S2cell_250_B"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "fruitfly_WT"
save_prefix <- "fruitfly_WT"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "fruitfly_luo2018_S2_normal"
save_prefix <- "fruitfly_luo2018_S2_normal"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "fruitfly_luo2018_S2_DMSO"
save_prefix <- "fruitfly_luo2018_S2_DMSO"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "fruitfly_SRR3031135"
save_prefix <- "fruitfly_SRR3031135"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

################################################################################################################

##### ZEBRAFISH #####
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.81_chr.gtf")), format="gtf")
#txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.84.chr.gtf")), format="gtf")
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

# make empty granges
starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                    IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                    Rle(txLengths$strand, txLengths$tx_len))

empty_gr$riboseq <- rep(0, length(empty_gr))

##########
wig_path <- "/Volumes/USELESS/STALLING/wigs/zebrafish"
pattern <- "24cell"
save_prefix <- "zebrafish_24cell"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "256cell"
save_prefix <- "zebrafish_256cell"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "Dome_trimmed"
save_prefix <- "zebrafish_Dome"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_2h"
save_prefix <- "zebrafish_2h"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_subtelny_2h"
save_prefix <- "zebrafish_subtelny_2h"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_subtelny_4h"
save_prefix <- "zebrafish_subtelny_4h"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_subtelny_6h"
save_prefix <- "zebrafish_subtelny_6h"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_beaudoin_2h"
save_prefix <- "zebrafish_beaudoin_2h"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_2h_1"
save_prefix <- "zebrafish_bazzini2012_2h_1"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_2h_2"
save_prefix <- "zebrafish_bazzini2012_2h_2"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_4h_1"
save_prefix <- "zebrafish_bazzini2012_4h_1"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_4h_2"
save_prefix <- "zebrafish_bazzini2012_4h_2"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_6h_1"
save_prefix <- "zebrafish_bazzini2012_6h_1"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2012_6h_2"
save_prefix <- "zebrafish_bazzini2012_6h_2"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_dome"
save_prefix <- "zebrafish_dome"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_shield"
save_prefix <- "zebrafish_shield"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "zebrafish_bazzini2014"
save_prefix <- "zebrafish_bazzini2014"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

################################################################################################################

##### MOUSE #####
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Mus_musculus.GRCm38.79.chr.gtf")), format="gtf")
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

# make empty granges
starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                    IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                    Rle(txLengths$strand, txLengths$tx_len))

empty_gr$riboseq <- rep(0, length(empty_gr))

##########
wig_path <- "/Volumes/USELESS/STALLING/wigs/mouse"
pattern <- "mouse_chx"
save_prefix <- "mouse_chx"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "mouse_emet"
save_prefix <- "mouse_emet"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "mouse_none"
save_prefix <- "mouse_none"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)

pattern <- "mouse_3T3"
save_prefix <- "mouse_3T3"
save_gr(empty_gr, pattern, wig_path, exons, save_prefix)
################################################################################################################

##### HUMAN #####
txdb_human <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")



# wig_path <- "/Volumes/USELESS/STALLING/wigs"
# import.wig("/Volumes/USELESS/STALLING/wigs/yeast/gb15_meio1-21-forward.wig")

###############################################################################################################
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
gr_path <- "/Volumes/USELESS/STALLING/gr"
wig_path <- "/Volumes/USELESS/STALLING/wigs"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  exons <- exonsBy(txdb, by="tx", use.names=TRUE)
  
  # make empty granges
  starts <- unlist(sapply(txLengths$tx_len, function(x){1:x}))
  empty_gr <- GRanges(seqnames=Rle(txLengths$tx_name, txLengths$tx_len),
                      IRanges(starts, width=rep(1, sum(txLengths$tx_len))),
                      Rle(txLengths$strand, txLengths$tx_len))
  
  empty_gr$riboseq <- rep(0, length(empty_gr))
  
  libs <- list.files(path = file.path(wig_path, org), pattern = paste0("^", org)) ## get org-specific libraries
  patterns <- unique(sapply(libs, function(x)substr(x, 1, (nchar(x)-15))))
  for (j in 1:length(patterns)) {
    pattern <- patterns[j]
    save_prefix <- pattern
    save_gr(empty_gr, pattern, file.path(wig_path, org), exons, save_prefix)
  }
  
}




