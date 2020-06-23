### stall sites on GRCz10 (for gollum)

library(GenomicFeatures)
library(rtracklayer)

path_bw <- "/Volumes/USELESS/OUT/zebrafish_our_out/single/from_gollum/shoelaces_GRCz10/bigwig"

cell24_fwd = file.path(path_bw, "cell24_fwd.bw")
cell24_rev = file.path(path_bw, "cell24_rev.bw")
cell256_fwd = file.path(path_bw, "cell256_fwd.bw")
cell256_rev = file.path(path_bw, "cell256_rev.bw")
cell1K_fwd = file.path(path_bw, "cell1K_fwd.bw")
cell1K_rev = file.path(path_bw, "cell1K_rev.bw")
dome_fwd = file.path(path_bw, "dome_fwd.bw")
dome_rev = file.path(path_bw, "dome_rev.bw")
shield_fwd = file.path(path_bw, "shield_fwd.bw")
shield_rev = file.path(path_bw, "shield_rev.bw")
bud_fwd = file.path(path_bw, "bud_fwd.bw")
bud_rev = file.path(path_bw, "bud_rev.bw")
hpf28_fwd = file.path(path_bw, "hpf28_fwd.bw")
hpf28_rev = file.path(path_bw, "hpf28_rev.bw")
dpf5_fwd = file.path(path_bw, "dpf5_fwd.bw")
dpf5_rev = file.path(path_bw, "dpf5_rev.bw")

bw_fw <- import.bw(cell24_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(cell24_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
cell24 <- c(bw_fw, bw_rv)

bw_fw <- import.bw(cell256_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(cell256_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
cell256 <- c(bw_fw, bw_rv)

bw_fw <- import.bw(cell1K_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(cell1K_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
cell1K <- c(bw_fw, bw_rv)

bw_fw <- import.bw(dome_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(dome_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
dome <- c(bw_fw, bw_rv)

bw_fw <- import.bw(shield_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(shield_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
shield <- c(bw_fw, bw_rv)

bw_fw <- import.bw(bud_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(bud_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
bud <- c(bw_fw, bw_rv)

bw_fw <- import.bw(hpf28_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(hpf28_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
hpf28 <- c(bw_fw, bw_rv)

bw_fw <- import.bw(dpf5_fwd)
strand(bw_fw) <- rep("+", length(bw_fw))
bw_rv <- import.bw(dpf5_rev)
strand(bw_rv) <- rep("-", length(bw_rv))
dpf5 <- c(bw_fw, bw_rv)

# txDb
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.GRCz10.81_chr.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
#u <- unlist(cds)
#u <- GRanges(seqnames = names(u), ranges(u), strand = strand(u))
#cds <- split(u, seqnames(u))
cds <- cds[order(names(cds))]


txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]




############
# sub CDS
LIB <- cell24
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
cell24_zscores <- zscores_he

# sub CDS
LIB <- cell256
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
cell256_zscores <- zscores_he

# sub CDS
LIB <- cell1K
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
cell1K_zscores <- zscores_he

# sub CDS
LIB <- dome
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
dome_zscores <- zscores_he

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

# sub CDS
LIB <- bud
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
bud_zscores <- zscores_he

# sub CDS
LIB <- hpf28
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
hpf28_zscores <- zscores_he

# sub CDS
LIB <- dpf5
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
dpf5_zscores <- zscores_he


####### get peaks
peaks_cell24 <- cell24_zscores[sapply(cell24_zscores, function(x){sum(x > 8)}) > 0]
peaks_cell24 <- sapply(peaks_cell24, function(x){which(x > 8)})

peaks_cell256 <- cell256_zscores[sapply(cell256_zscores, function(x){sum(x > 8)}) > 0]
peaks_cell256 <- sapply(peaks_cell256, function(x){which(x > 8)})

peaks_cell1K <- cell1K_zscores[sapply(cell1K_zscores, function(x){sum(x > 8)}) > 0]
peaks_cell1K <- sapply(peaks_cell1K, function(x){which(x > 8)})

peaks_dome <- dome_zscores[sapply(dome_zscores, function(x){sum(x > 8)}) > 0]
peaks_dome <- sapply(peaks_dome, function(x){which(x > 8)})

peaks_shield <- shield_zscores[sapply(shield_zscores, function(x){sum(x > 8)}) > 0]
peaks_shield <- sapply(peaks_shield, function(x){which(x > 8)})

peaks_bud <- bud_zscores[sapply(bud_zscores, function(x){sum(x > 8)}) > 0]
peaks_bud <- sapply(peaks_bud, function(x){which(x > 8)})

peaks_hpf28 <- hpf28_zscores[sapply(hpf28_zscores, function(x){sum(x > 8)}) > 0]
peaks_hpf28 <- sapply(peaks_hpf28, function(x){which(x > 8)})

peaks_dpf5 <- dpf5_zscores[sapply(dpf5_zscores, function(x){sum(x > 8)}) > 0]
peaks_dpf5 <- sapply(peaks_dpf5, function(x){which(x > 8)})


peaks <- list(cell24=peaks_cell24, cell256=peaks_cell256, cell1K=peaks_cell1K, dome=peaks_dome, shield=peaks_shield, bud=peaks_bud, hpf28=peaks_hpf28, dpf5=peaks_dpf5)

save(peaks, file="/Users/kasia/Documents/PhD/scripts/STALLING/stall_sites_GRCz10/peaks_GRCz10.Rsave")

