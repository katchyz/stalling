### for SNP analysis:
library(rtracklayer)

# granges with stall site codons
h_tx_ss <- css[!is.na(css$h_ss),][,10:11]
h_tx_ss$end <- h_tx_ss$h_ss + 3
# h_tx_ss$h_ss <- h_tx_ss$h_ss + 1
# minus <- ss_genomic[strand(ss_genomic) == "-"]$seqnames
# h_tx_ss[h_tx_ss$h_tx %in% minus,]$h_ss <- h_tx_ss[h_tx_ss$h_tx %in% minus,]$h_ss + 1
# h_tx_ss[h_tx_ss$h_tx %in% minus,]$end <- h_tx_ss[h_tx_ss$h_tx %in% minus,]$end + 1

cds_h <- cdsBy(txdb_h, by="tx", use.names=TRUE)
cds <- cds_h[names(cds_h) %in% h_tx_ss$h_tx] # 1729 tx

gr <- makeGRangesFromDataFrame(h_tx_ss,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="h_tx",
                               start.field="h_ss",
                               end.field="end",
                               starts.in.df.are.0based=TRUE)

ss_genomic <- mapFromTranscripts(gr, cds)
ss_genomic$seqnames <- as.character(seqnames(gr))
export.bed(ss_genomic, "/Volumes/USELESS/STALLING/SNP/human.bed")

# bam file(s): check peak files whether they occur in a given library

h_tx_ss$h_ss <- h_tx_ss$h_ss + 1
h_tx_ss$tx_ss <- paste0(h_tx_ss$h_tx, "_", h_tx_ss$h_ss)

# h_tx_ss <- h_tx_ss[order(h_tx_ss$h_tx),]

elementMetadata(ss_genomic)$tx <- h_tx_ss$tx_ss

save(ss_genomic, file = "/Volumes/USELESS/STALLING/SNP/human_genomic_granges.Rsave")

peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median/"
libs <- list.files(path = peaks_path, pattern = "human")

h_tx_ss$h1 <- rep(0, nrow(h_tx_ss))
h_tx_ss$h2 <- rep(0, nrow(h_tx_ss))
h_tx_ss$h3 <- rep(0, nrow(h_tx_ss))
h_tx_ss$h4 <- rep(0, nrow(h_tx_ss))

for (i in 1:length(libs)) {
  ## load file
  load(file = c(file.path(peaks_path, libs[i])))
  peaks$utr5_len <- sapply(peaks$seqnames, function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
  peaks$ss <- peaks$start - peaks$utr5_len
  
  peaks$tx_ss <- paste0(peaks$seqnames, "_", peaks$ss)
  print(libs[i])
  print(sum(h_tx_ss$tx_ss %in% peaks$tx_ss))
  
  h_tx_ss[h_tx_ss$tx_ss %in% peaks$tx_ss,][,(i+3)] <- 1
  
}

save(h_tx_ss, file = "/Volumes/USELESS/STALLING/SNP/human_table_bam.Rsave")


######## control SNPs

h_tx_ss <- css[!is.na(css$h_ss),][,10:11]
h_tx_ss$h_ss <- h_tx_ss$h_ss + 1
h_tx_ss$end <- h_tx_ss$h_ss + 2

h_tx_ss$cds_len <- sapply(as.character(h_tx_ss$h_tx), function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

gr_ss <- makeGRangesFromDataFrame(h_tx_ss,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="h_tx",
                               start.field="h_ss",
                               end.field="end",
                               starts.in.df.are.0based=FALSE)

h_tx_ss$random_start <- sapply(h_tx_ss$cds_len, function(x){sample(seq(16,(x-5),3),1)})
h_tx_ss[abs(h_tx_ss$h_ss - h_tx_ss$random_start) < 4,]$random_start <- 
  sapply(h_tx_ss[abs(h_tx_ss$h_ss - h_tx_ss$random_start) < 4,]$cds_len, function(x){sample(seq(16,x,3),1)})

h_tx_ss$random_end <- h_tx_ss$random_start + 2

gr_random <- makeGRangesFromDataFrame(h_tx_ss,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="h_tx",
                                  start.field="random_start",
                                  end.field="random_end",
                                  starts.in.df.are.0based=FALSE)


ss_genomic <- mapFromTranscripts(gr_ss, cds)
ss_genomic$seqnames <- as.character(seqnames(gr_ss))
export.bed(ss_genomic, "/Volumes/USELESS/STALLING/SNP/human_CSS_SNP.bed")

random_genomic <- mapFromTranscripts(gr_random, cds)
random_genomic$seqnames <- as.character(seqnames(gr_random))
export.bed(random_genomic, "/Volumes/USELESS/STALLING/SNP/human_random_SNP.bed")


