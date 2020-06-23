### peak stats
# how many in each library
# how many repeated, between two given libs fro an organism
# decide on final libs, different experiments?
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(plyr)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

# get the longest tx per gene only

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

#summary table
summary_table <- data.table(library1=c("lib1"), library2=c("lib2"), lib1=c(0), lib2=c(0), common=c(0))

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
  #dt <- data.table(seqnames=c("transcript"), start=c(0), library=c("library"))
  ### longest transcript per gene
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  ##
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  ## get all possible couples of libraries
  couples <- combn(libs, 2)
  for (j in 1:ncol(couples)) {
    peaks1 <- get(load(file = c(file.path(peaks_path, couples[1,j]))))
    peaks2 <- get(load(file = c(file.path(peaks_path, couples[2,j]))))
    # get longest transcript only
    peaks1 <- peaks1[peaks1$seqnames %in% rownames(txlen)]
    peaks2 <- peaks2[peaks2$seqnames %in% rownames(txlen)]
    peaks1 <- peaks1[,1:2]
    peaks2 <- peaks2[,1:2]
    # merge
    p1 <- paste(peaks1$seqnames, peaks1$start, sep="_")
    p2 <- paste(peaks2$seqnames, peaks2$start, sep="_")
    # summary
    st <- data.table(library1=couples[1,j], library2=couples[2,j], lib1=nrow(peaks1), lib2=nrow(peaks2),
                     common = sum(p1 %in% p2))
    summary_table <- rbind(summary_table, st)
  }
}
summary_table <- summary_table[-1,] # delete first row (was just a placeholder)
save(summary_table, file = c(file.path(peaks_path, "summary_table.Rsave")))

###########
# load peaks, get only longest transcript, redo consensus
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  nlib <- length(list.files(path = peaks_path, pattern = paste0("^", org)))
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  setDT(consensus)
  consensus <- consensus[consensus$start > 15]
  ### longest tx
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  ##
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  ###
  consensus <- consensus[consensus$seqnames %in% rownames(txlen)]
  sum(consensus$freq > 1)
  sum(consensus$freq == nlib)
  #consensus <- consensus[consensus$freq == nlib,]
}


