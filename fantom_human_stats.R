### reads mapping to CDS, UTRs
library(GenomicFeatures)
library(GenomicAlignments)
library(plyr)

# load gtf, get cds, utr5
txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# get transcript with longest cds per gene
txlen <- arrange(txLengths, gene_id, desc(cds_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]
# features
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)

cds <- cds[names(cds) %in% rownames(txlen)]
utr5 <- utr5[names(utr5) %in% rownames(txlen)]
utr3 <- utr3[names(utr3) %in% rownames(txlen)]

f = "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/stats.txt"
cat(c("study", "total", "CDS", "utr5", "utr3", "\n"), file=f, sep="\t")

# load bam
bam_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/human"
libs <- list.files(path = bam_path, pattern = "\\.bam$")

for (i in 1:length(libs)) {
  ## load file
  file = c(file.path(bam_path, libs[i]))
  alignments <- readGAlignments(file)
  aln_cds <- subsetByOverlaps(alignments, cds)
  aln_utr5 <- subsetByOverlaps(alignments, utr5)
  aln_utr3 <- subsetByOverlaps(alignments, utr3)
  cat(c(libs[i], length(alignments), length(aln_cds), length(aln_utr5), length(aln_utr3), "\n"), file=f, append=TRUE)
  
  flengths <- unique(qwidth(alignments))
  for (j in 1:length(flengths)) {
    laln <- length(alignments[qwidth(alignments) == flengths[j]])
    laln_cds <- length(aln_cds[qwidth(aln_cds) == flengths[j]])
    laln_utr5 <- length(aln_utr5[qwidth(aln_utr5) == flengths[j]])
    laln_utr3 <- length(aln_utr3[qwidth(aln_utr3) == flengths[j]])
    cat(c(flengths[j], laln, laln_cds, laln_utr5, laln_utr3, "\n"), file=f, append=TRUE)
  }
}



