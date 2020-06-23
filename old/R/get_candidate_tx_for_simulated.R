### get good candidate transcripts for simulated bam

## read ribo-seq library
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(GeneCycle)

path_bw <- "/Volumes/USELESS/OUT/zebrafish_our_out/single/from_gollum/shoelaces_GRCz10/bigwig"

shield_fwd = file.path(path_bw, "shield_fwd.bw")
shield_rev = file.path(path_bw, "shield_rev.bw")

path_bw <- "/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/by_length/bw"

shield_fwd = file.path(path_bw, "Shield_trimmed-28-forward.bw")
shield_rev = file.path(path_bw, "Shield_trimmed-28-reverse.bw")

bw_fw <- import.bw(shield_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(shield_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
shield <- c(bw_fw, bw_rv)

# txDb
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.GRCz10.81_chr.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# sub CDS
LIB <- shield
LIB_cds <- mapToTranscripts(LIB, cds, ignore.strand = FALSE)
mcols(LIB_cds) <- cbind(mcols(LIB_cds), DataFrame(LIB[LIB_cds$xHits]))
LIB_cds <- split(LIB_cds, seqnames(LIB_cds))
LIB_cds <- LIB_cds[order(names(LIB_cds))]
# fill vectors with 0s
seqlengths(LIB_cds)[order(names(seqlengths(LIB_cds)))] <- 
  txLengths[rownames(txLengths) %in% names(LIB_cds),]$cds_len
LIB_cov_cds <- coverage(LIB_cds, weight=unlist(LIB_cds)$score)
LIB_cov_cds <- sapply(LIB_cov_cds, function(x){as.vector(x)})
LIB_cov_cds <- LIB_cov_cds[sapply(LIB_cov_cds, function(x){sum(x) > 10})]
LIB_codons <- sapply(LIB_cov_cds, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
he <- LIB_codons[sapply(LIB_codons, function(x){median(x) > 0})]

zscores_he <- lapply(he, function(x){c(0,as.vector(scale(x[2:(length(x)-1)])),0)})
shield_zscores <- zscores_he

#############
## get exons
## tx where highest peaks are around start and stop

LIB <- shield
LIB_exons <- mapToTranscripts(LIB, exons, ignore.strand = FALSE)
mcols(LIB_exons) <- cbind(mcols(LIB_exons), DataFrame(LIB[LIB_exons$xHits]))
LIB_exons <- split(LIB_exons, seqnames(LIB_exons))
LIB_exons <- LIB_exons[order(names(LIB_exons))]
# fill vectors with 0s
seqlengths(LIB_cds)[order(names(seqlengths(LIB_cds)))] <- 
  txLengths[rownames(txLengths) %in% names(LIB_cds),]$tx_len
LIB_cov_exons <- coverage(LIB_exons, weight=unlist(LIB_exons)$score)
LIB_cov_exons <- sapply(LIB_cov_exons, function(x){as.vector(x)})
LIB_cov_exons <- LIB_cov_exons[sapply(LIB_cov_exons, function(x){sum(x) > 100})]
LIB_cov_exons <- LIB_cov_exons[order(names(LIB_cov_exons))]

### are the highest peaks around start/stop?
utr5_len <- txLengths[names(LIB_cov_exons),]$utr5_len
cds_len <- txLengths[names(LIB_cov_exons),]$cds_len

thr <- sapply(LIB_cov_exons, function(x){head(sort(x, decreasing = TRUE), 3)[3]})
highest_peaks <- mapply(function(x,y){which(x > y)}, LIB_cov_exons, thr)

start <- sapply(utr5_len, function(x){seq(x-2,x+3)})
stop <- sapply(utr5_len+cds_len, function(x){seq(x-5,x+3)})
start_stop <- rbind(start, stop)

start <- as.list(data.frame(start))
stop <- as.list(data.frame(stop))
start_stop <- as.list(data.frame(start_stop))

sum(mapply(function(x,y){sum(x %in% y) > 1}, highest_peaks, start_stop))

hp_start <- names(highest_peaks)[mapply(function(x,y){sum(x %in% y) > 0}, highest_peaks, start)]
hp_stop <- names(highest_peaks)[mapply(function(x,y){sum(x %in% y) > 0}, highest_peaks, stop)]

h <- hp_start[hp_start %in% hp_stop]


# check for periodicity
meta <- LIB_cov_cds[[h[16]]][1:150]

amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell256, TC.control") + xlim(0, 8)


tx_start <- c()
meta_cum <- rep(0,150)
for (i in 1:length(hp_start)) {
  if (length(LIB_cov_cds[[hp_start[i]]]) > 150) {
    meta <- LIB_cov_cds[[hp_start[i]]][1:150]
    meta_cum <- meta_cum + meta
  }
  if (length(meta) == 150) {
    amplitudes <- abs(fft(meta))
    amp <- amplitudes[2:(length(amplitudes)/2+1)]
    periods <- 1/periodogram(meta, method = "clone")$freq
    if (periods[which.max(amp)] == 3) {
      print(hp_start[i])
      tx_start <- c(tx_start, hp_start[i])
    }
  }
}


tx_stop <- c()
meta_cum <- rep(0,150)
for (i in 1:length(hp_stop)) {
  if (length(LIB_cov_cds[[hp_stop[i]]]) > 150) {
    meta <- LIB_cov_cds[[hp_stop[i]]][1:150]
    meta_cum <- meta_cum + meta
  }
  if (length(meta) == 150) {
    amplitudes <- abs(fft(meta))
    amp <- amplitudes[2:(length(amplitudes)/2+1)]
    periods <- 1/periodogram(meta, method = "clone")$freq
    if (periods[which.max(amp)] == 3) {
      print(hp_stop[i])
      tx_stop <- c(tx_stop, hp_stop[i])
    }
  }
}


gene_start <- txLengths[tx_start,]$gene_id
gene_stop <- txLengths[tx_stop,]$gene_id

for (i in 1:length(gene_start)) {
  print(gene_start[i])
}

for (i in 1:length(gene_stop)) {
  print(gene_stop[i])
}


amplitudes <- abs(fft(meta_cum))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta_cum, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell256, TC.control") + xlim(0, 8)



