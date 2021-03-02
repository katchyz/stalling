# k-mers for Carl

library(data.table)
library(GenomicFeatures)
library(seqinr)

# stall sites
css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
human_ss <- data.table(tx = css[!is.na(css$h_ss),]$h_tx, gene = css[!is.na(css$h_ss),]$gene,
                       nt_cds = css[!is.na(css$h_ss),]$h_ss + 1)

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]

human_ss$cds_len <- sapply(human_ss$tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

# control positions
### CONTROL
utx <- unique(as.character(human_ss$tx))
random_pool <- list()
for (i in 1:length(utx)) {
  txn <- utx[i]
  cds_len_nt <- human_ss[human_ss$tx == txn,]$cds_len[1]
  around_ss <- c()
  for (ss in human_ss[as.character(human_ss$tx) == txn,]$nt_cds) {
    around_ss <- c(around_ss, c((ss-99):(ss+101)))
  }
  pool <- c(101:(cds_len_nt-101))[!c(101:(cds_len_nt-101)) %in% around_ss]
  random_pool[[txn]] <- pool
}

random_pool <- random_pool[sapply(random_pool, function(x){length(x)}) > 0]
human_ss_random <- human_ss[human_ss$tx %in% names(random_pool)]

human_ss_random$random1 <- sapply(as.character(human_ss_random$tx), function(x){sample(random_pool[[x]],1)})
human_ss_random$random2 <- sapply(as.character(human_ss_random$tx), function(x){sample(random_pool[[x]],1)})
human_ss_random <- human_ss_random[human_ss_random$nt_cds > 100]
human_ss_random <- human_ss_random[human_ss_random$random1 > 100]
# human_ss_random <- human_ss_random[human_ss_random$cds_len - human_ss_random$random1 > 101]
# human_ss_random <- human_ss_random[!human_ss_random$tx == "ENST00000496455"]

# fasta
fasta_cds <- read.fasta("../../DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz")

css_neigh <- t(mapply(function(x,y){fasta_cds[[x]][(y-99):(y+101)]}, as.character(human_ss_random$tx),
                      human_ss_random$nt_cds))
tn <- rownames(css_neigh)
css_neigh <- as.data.table(css_neigh)

css_random <- t(mapply(function(x,y){fasta_cds[[x]][(y-99):(y+101)]}, as.character(human_ss_random$tx),
                       human_ss_random$random1))
tr <- rownames(css_random)
css_random <- as.data.table(css_random)

css201 <- data.table(tx = tn, seq = do.call(paste0, css_neigh[,2:201])) # neighbourhood of stall sites
rss201 <- data.table(tx = tr, seq = do.call(paste0, css_random[,2:201])) # neighbourhood of control positions

kmer <- data.table(tx = css201$tx, ss = human_ss_random$nt_cds, css_seq = css201$seq, control_seq = rss201$seq)
# save(kmer, file = "../../sequence/kmer.Rsave")

#######################


