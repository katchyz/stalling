# check riboseq before and after for Human (H1-H4)
library(GenomicFeatures)
library(ggplot2)

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]

ht <- get(load(file = "../../DATA/termination/human_table_bam.Rsave"))
h2 <- ht[ht$h2 == 1,]
h2$utr5_len <- sapply(h2$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
h2$cds_len <- sapply(h2$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
h2$start_cdna <- h2$h_ss + h2$utr5_len
h2$dist_to_end <- h2$cds_len - h2$h_ss - 2

h100 <- h2[h2$h_ss > 145 & h2$dist_to_end > 106,]

load(file = "../../DATA/gr/human_stern_NODRUG.Rsave") # gr
gr <- split(gr, seqnames(gr))

gr <- gr[names(gr) %in% as.character(h2$h_tx)]

dt_log2 <- data.table(tx = "tx", start = 0, log2 = 0)
for (i in 1:nrow(h100)) {
  tx <- h100[i,]$h_tx
  start <- h100[i,]$start_cdna
  
  before <- gr[[tx]]$riboseq[(start-100):(start-1)]
  after <- gr[[tx]]$riboseq[(start+3):(start+102)]
  
  if (sum(after) > 0) {
    log2ratio <- log2(sum(before)/sum(after))
  } else {
    log2ration <- NA
  }
  
  dt <- data.table(tx = tx, start = start, log2 = log2ratio)
  dt_log2 <- rbind(dt_log2, dt)
  
}

save(dt_log2, file = "../../DATA/termination/log2ratio_h2.Rsave")
ggplot(dt_log2, aes(x = log2)) + geom_histogram() + theme_classic()

### whole before and after
h_whole <- h2[h2$h_ss > 90 & h2$dist_to_end > 45,]

dt_log2 <- data.table(tx = "tx", start = 0, log2 = 0)
for (i in 1:nrow(h_whole)) {
  tx <- h_whole[i,]$h_tx
  start <- h_whole[i,]$start_cdna
  utr5 <- h_whole[i,]$utr5_len
  to_end <- h_whole[i,]$dist_to_end
  
  before <- gr[[tx]]$riboseq[(start-utr5-45):(start-1)]
  after <- gr[[tx]]$riboseq[(start+3):(start+to_end)]
  
  if (sum(after) > 0) {
    log2ratio <- log2(mean(before)/mean(after))
  } else {
    log2ratio <- NA
  }
  
  dt <- data.table(tx = tx, start = start, log2 = log2ratio)
  dt_log2 <- rbind(dt_log2, dt)
}

save(dt_log2, file = "../../DATA/termination/log2ratio_h2_whole.Rsave")
ggplot(dt_log2, aes(x = log2)) + geom_histogram() + theme_classic()
# ggsave("log2ratio_h2_whole.png")

### H1
load(file = "../../DATA/gr/human_fibroblasts_stern2012.Rsave") # h1
gr <- split(gr, seqnames(gr))
h1 <- ht[ht$h1 == 1,]
h1$utr5_len <- sapply(h1$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
h1$cds_len <- sapply(h1$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
h1$start_cdna <- h1$h_ss + h1$utr5_len
h1$dist_to_end <- h1$cds_len - h1$h_ss - 2
h_whole <- h1[h1$h_ss > 90 & h1$dist_to_end > 45,]

gr <- gr[names(gr) %in% as.character(h1$h_tx)]

dt_log2 <- data.table(tx = "tx", start = 0, log2 = 0)
for (i in 1:nrow(h_whole)) {
  tx <- h_whole[i,]$h_tx
  start <- h_whole[i,]$start_cdna
  utr5 <- h_whole[i,]$utr5_len
  to_end <- h_whole[i,]$dist_to_end
  
  before <- gr[[tx]]$riboseq[(start-utr5-45):(start-1)]
  after <- gr[[tx]]$riboseq[(start+3):(start+to_end)]
  
  if (sum(after) > 0) {
    log2ratio <- log2(mean(before)/mean(after))
  } else {
    log2ratio <- NA
  }
  
  dt <- data.table(tx = tx, start = start, log2 = log2ratio)
  dt_log2 <- rbind(dt_log2, dt)
}
save(dt_log2, file = "../../DATA/termination/log2ratio_h1_whole.Rsave")
ggplot(dt_log2, aes(x = log2)) + geom_histogram() + theme_classic()
# ggsave("log2ratio_h1_whole.png")


### H3
load(file = "../../DATA/gr/human_HeLa_stumpf2013.Rsave") # h3
gr <- split(gr, seqnames(gr))
h3 <- ht[ht$h3 == 1,]
h3$utr5_len <- sapply(h3$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
h3$cds_len <- sapply(h3$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
h3$start_cdna <- h3$h_ss + h3$utr5_len
h3$dist_to_end <- h3$cds_len - h3$h_ss - 2
h_whole <- h3[h3$h_ss > 90 & h3$dist_to_end > 45,]

gr <- gr[names(gr) %in% as.character(h3$h_tx)]

dt_log2 <- data.table(tx = "tx", start = 0, log2 = 0)
for (i in 1:nrow(h_whole)) {
  tx <- h_whole[i,]$h_tx
  start <- h_whole[i,]$start_cdna
  utr5 <- h_whole[i,]$utr5_len
  to_end <- h_whole[i,]$dist_to_end
  
  before <- gr[[tx]]$riboseq[(start-utr5-45):(start-1)]
  after <- gr[[tx]]$riboseq[(start+3):(start+to_end)]
  
  if (sum(after) > 0) {
    log2ratio <- log2(mean(before)/mean(after))
  } else {
    log2ratio <- NA
  }
  
  dt <- data.table(tx = tx, start = start, log2 = log2ratio)
  dt_log2 <- rbind(dt_log2, dt)
}
save(dt_log2, file = "../../DATA/termination/log2ratio_h3_whole.Rsave")
ggplot(dt_log2, aes(x = log2)) + geom_histogram() + theme_classic()
# ggsave("log2ratio_h3_whole.png")



### H4
load(file = "../../DATA/gr/human_HEK_subtelny2014.Rsave") # h4
gr <- split(gr, seqnames(gr))
h4 <- ht[ht$h4 == 1,]
h4$utr5_len <- sapply(h4$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
h4$cds_len <- sapply(h4$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
h4$start_cdna <- h4$h_ss + h4$utr5_len
h4$dist_to_end <- h4$cds_len - h4$h_ss - 2
h_whole <- h4[h4$h_ss > 90 & h4$dist_to_end > 45,]

gr <- gr[names(gr) %in% as.character(h4$h_tx)]

dt_log2 <- data.table(tx = "tx", start = 0, log2 = 0)
for (i in 1:nrow(h_whole)) {
  tx <- h_whole[i,]$h_tx
  start <- h_whole[i,]$start_cdna
  utr5 <- h_whole[i,]$utr5_len
  to_end <- h_whole[i,]$dist_to_end
  
  before <- gr[[tx]]$riboseq[(start-utr5-45):(start-1)]
  after <- gr[[tx]]$riboseq[(start+3):(start+to_end)]
  
  if (sum(after) > 0) {
    log2ratio <- log2(mean(before)/mean(after))
  } else {
    log2ratio <- NA
  }
  
  dt <- data.table(tx = tx, start = start, log2 = log2ratio)
  dt_log2 <- rbind(dt_log2, dt)
}
save(dt_log2, file = "../../DATA/termination/log2ratio_h4_whole.Rsave")
ggplot(dt_log2, aes(x = log2)) + geom_histogram() + theme_classic()
# ggsave("log2ratio_h4_whole.png")


###################################################################################
# load data tables
# check the tx with large change
# check in homologs in other organisms, if it looks promising

h1 <- get(load(file = "../../DATA/termination/log2ratio_h1_whole.Rsave"))
h2 <- get(load(file = "../../DATA/termination/log2ratio_h2_whole.Rsave"))
h3 <- get(load(file = "../../DATA/termination/log2ratio_h3_whole.Rsave"))
h4 <- get(load(file = "../../DATA/termination/log2ratio_h4_whole.Rsave"))

css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})

css[css$h_tx %in% h1[h1$log2 > 2]$tx,]
css[css$h_tx %in% h2[h2$log2 > 2]$tx,]
css[css$h_tx %in% h3[h3$log2 > 2]$tx,]
css[css$h_tx %in% h4[h4$log2 > 2]$tx,]

ggplot(h1, aes(x = log2)) + geom_histogram() + theme_classic() + geom_vline(xintercept=2, linetype="dashed", color = "red")
# ggsave("red_h1.png")
ggplot(h2, aes(x = log2)) + geom_histogram() + theme_classic() + geom_vline(xintercept=2, linetype="dashed", color = "red")
# ggsave("red_h2.png")
ggplot(h3, aes(x = log2)) + geom_histogram() + theme_classic() + geom_vline(xintercept=2, linetype="dashed", color = "red")
# ggsave("red_h3.png")
ggplot(h4, aes(x = log2)) + geom_histogram() + theme_classic() + geom_vline(xintercept=2, linetype="dashed", color = "red")
# ggsave("red_h4.png")




