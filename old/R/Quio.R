### stalling (in genes from Quio)
library(GenomicFeatures)
library(rtracklayer)

## get ribo cov (on Zv9)
path_bw <- "/Volumes/USELESS/OUT/zebrafish_our_out/single"

cell24_fwd = file.path(path_bw, "cell24_fwd.bw")
cell24_rev = file.path(path_bw, "cell24_rev.bw")
cell256_fwd = file.path(path_bw, "cell256_fwd.bw")
cell256_rev = file.path(path_bw, "cell256_rev.bw")
cell1K_fwd = file.path(path_bw, "cell1K_fwd.bw")
cell1K_rev = file.path(path_bw, "cell1K_rev.bw")

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

# txDb
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

################ cell24
# sub CDS
cell24_cds <- mapToTranscripts(cell24, cds, ignore.strand = FALSE)
mcols(cell24_cds) <- cbind(mcols(cell24_cds), DataFrame(cell24[cell24_cds$xHits]))
cell24_cds <- split(cell24_cds, seqnames(cell24_cds))
cell24_cds <- cell24_cds[order(names(cell24_cds))]
# sub exons
cell24_exons <- mapToTranscripts(cell24, exons, ignore.strand = FALSE)
mcols(cell24_exons) <- cbind(mcols(cell24_exons), DataFrame(cell24[cell24_exons$xHits]))
cell24_exons <- split(cell24_exons, seqnames(cell24_exons))
cell24_exons <- cell24_exons[order(names(cell24_exons))]

# fill vectors with 0s
seqlengths(cell24_cds)[order(names(seqlengths(cell24_cds)))] <- 
  txLengths[rownames(txLengths) %in% names(cell24_cds),]$cds_len

seqlengths(cell24_exons)[order(names(seqlengths(cell24_exons)))] <- 
  txLengths[rownames(txLengths) %in% names(cell24_exons),]$tx_len

#coverage(ribo_cds, weight=unlist(ribo_cds)$score)$ENSDART00000000004
cell24_cov_cds <- coverage(cell24_cds, weight=unlist(cell24_cds)$score)
cell24_cov_exons <- coverage(cell24_exons, weight=unlist(cell24_exons)$score)

################ cell256
# sub CDS
cell256_cds <- mapToTranscripts(cell256, cds, ignore.strand = FALSE)
mcols(cell256_cds) <- cbind(mcols(cell256_cds), DataFrame(cell256[cell256_cds$xHits]))
cell256_cds <- split(cell256_cds, seqnames(cell256_cds))
cell256_cds <- cell256_cds[order(names(cell256_cds))]
# sub exons
cell256_exons <- mapToTranscripts(cell256, exons, ignore.strand = FALSE)
mcols(cell256_exons) <- cbind(mcols(cell256_exons), DataFrame(cell256[cell256_exons$xHits]))
cell256_exons <- split(cell256_exons, seqnames(cell256_exons))
cell256_exons <- cell256_exons[order(names(cell256_exons))]

# fill vectors with 0s
seqlengths(cell256_cds)[order(names(seqlengths(cell256_cds)))] <- 
  txLengths[rownames(txLengths) %in% names(cell256_cds),]$cds_len

seqlengths(cell256_exons)[order(names(seqlengths(cell256_exons)))] <- 
  txLengths[rownames(txLengths) %in% names(cell256_exons),]$tx_len

#coverage(ribo_cds, weight=unlist(ribo_cds)$score)$ENSDART00000000004
cell256_cov_cds <- coverage(cell256_cds, weight=unlist(cell256_cds)$score)
cell256_cov_exons <- coverage(cell256_exons, weight=unlist(cell256_exons)$score)

################ cell1K
# sub CDS
cell1K_cds <- mapToTranscripts(cell1K, cds, ignore.strand = FALSE)
mcols(cell1K_cds) <- cbind(mcols(cell1K_cds), DataFrame(cell1K[cell1K_cds$xHits]))
cell1K_cds <- split(cell1K_cds, seqnames(cell1K_cds))
cell1K_cds <- cell1K_cds[order(names(cell1K_cds))]
# sub exons
cell1K_exons <- mapToTranscripts(cell1K, exons, ignore.strand = FALSE)
mcols(cell1K_exons) <- cbind(mcols(cell1K_exons), DataFrame(cell1K[cell1K_exons$xHits]))
cell1K_exons <- split(cell1K_exons, seqnames(cell1K_exons))
cell1K_exons <- cell1K_exons[order(names(cell1K_exons))]

# fill vectors with 0s
seqlengths(cell1K_cds)[order(names(seqlengths(cell1K_cds)))] <- 
  txLengths[rownames(txLengths) %in% names(cell1K_cds),]$cds_len

seqlengths(cell1K_exons)[order(names(seqlengths(cell1K_exons)))] <- 
  txLengths[rownames(txLengths) %in% names(cell1K_exons),]$tx_len

#coverage(ribo_cds, weight=unlist(ribo_cds)$score)$ENSDART00000000004
cell1K_cov_cds <- coverage(cell1K_cds, weight=unlist(cell1K_cds)$score)
cell1K_cov_exons <- coverage(cell1K_exons, weight=unlist(cell1K_exons)$score)

############## sum up cell24, cell256 and cell1K ..._cov_cds and _cov_exons
c24 <- cell24_cov_cds[names(cell24_cov_cds) %in% names(cell256_cov_cds)]
c256 <- cell256_cov_cds[names(cell256_cov_cds) %in% names(c24)]
cov_cds <- c24 + c256
cov_cds <- c(cov_cds, cell24_cov_cds[!(names(cell24_cov_cds) %in% names(cov_cds))])
cov_cds <- c(cov_cds, cell256_cov_cds[!(names(cell256_cov_cds) %in% names(cov_cds))])
c1K <- cell1K_cov_cds[names(cell1K_cov_cds) %in% names(cov_cds)]
cov_cds_c1K <- cov_cds[names(cov_cds) %in% names(c1K)]
cov_cds_c1K <- cov_cds_c1K[order(names(cov_cds_c1K))]
c1K <- c1K[order(names(c1K))]
coverage_cds <- cov_cds_c1K + c1K
coverage_cds <- c(coverage_cds, cov_cds[!(names(cov_cds) %in% names(coverage_cds))])
coverage_cds <- coverage_cds[order(names(coverage_cds))] #!!!!

c24 <- cell24_cov_exons[names(cell24_cov_exons) %in% names(cell256_cov_exons)]
c256 <- cell256_cov_exons[names(cell256_cov_exons) %in% names(c24)]
cov_exons <- c24 + c256
cov_exons <- c(cov_exons, cell24_cov_exons[!(names(cell24_cov_exons) %in% names(cov_exons))])
cov_exons <- c(cov_exons, cell256_cov_exons[!(names(cell256_cov_exons) %in% names(cov_exons))])
c1K <- cell1K_cov_exons[names(cell1K_cov_exons) %in% names(cov_exons)]
cov_exons_c1K <- cov_exons[names(cov_exons) %in% names(c1K)]
cov_exons_c1K <- cov_exons_c1K[order(names(cov_exons_c1K))]
c1K <- c1K[order(names(c1K))]
coverage_exons <- cov_exons_c1K + c1K
coverage_exons <- c(coverage_exons, cov_exons[!(names(cov_exons) %in% names(coverage_exons))])
coverage_exons <- coverage_exons[order(names(coverage_exons))] #!!!!

#############
## subset genes of interest (Degraged_uORF_genes.txt)
tx_gene <- read.table(file="/Users/kasia/Documents/PhD/scripts/STALLING/Quio/Degraded_uORF_genes.txt", header=T)

subset_cov_cds <- coverage_cds[names(coverage_cds) %in% tx_gene$Transcript]
subset_cov_cds <- subset_cov_cds[order(names(subset_cov_cds))]
subset_cov_exons <- coverage_exons[names(coverage_exons) %in% tx_gene$Transcript]
subset_cov_exons <- subset_cov_exons[order(names(subset_cov_exons))]


make_gr <- function(x,y){GRanges(seqnames = rep(y, length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep("*", length(x)),
                                   ribo_cov = as.vector(x))} # converting Rle into GRanges

gr_cds <- GRangesList(mapply(make_gr, subset_cov_cds, names(subset_cov_cds)))

gr_exons <- GRangesList(mapply(make_gr, subset_cov_exons, names(subset_cov_exons)))

## ?? check for uORFs ??
## calculate z-scores per 3nt (remove start and stop codon)

#sum every 3nt
ribo <- sapply(gr_cds, function(x){as.vector(tapply(x$ribo_cov[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})

## check coverage on CDS (median > 0)
sum(sapply(ribo, function(x){median(x) > 0}))

he <- ribo[sapply(ribo, function(x){median(x) > 0})] #highly expressed


#zscores_cds <- sapply(gr_cds, function(x){c(0,0,0,as.vector(scale(x$ribo_cov[4:(length(x)-4)])),0,0,0)})
zscores <- sapply(ribo, function(x){c(0,as.vector(scale(x[2:(length(x)-1)])),0)})

zscores_he <- sapply(he, function(x){c(0,as.vector(scale(x[2:(length(x)-1)])),0)})

## check for peaks (z-score threshold)
sapply(zscores, function(x){sum(x > 5)})
sapply(zscores_he, function(x){sum(x > 5)})


############# uORFs

tx_he <- gr_exons[names(gr_exons) %in% names(he)]



########### CONTROL on non-degraded transcripts
gr_cds_all <- GRangesList(mapply(make_gr, coverage_cds, names(coverage_cds)))
gr_exons_all <- GRangesList(mapply(make_gr, coverage_exons, names(coverage_exons)))

ribo_all <- sapply(gr_cds_all, function(x){as.vector(tapply(x$ribo_cov[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
## check coverage on CDS (median > 0)
sum(sapply(ribo_all, function(x){median(x) > 0}))
he_all <- ribo_all[sapply(ribo_all, function(x){median(x) > 0})] #highly expressed
#zscores <- sapply(ribo_all, function(x){c(0,as.vector(scale(x[2:(length(x)-1)])),0)})
zscores_he_all <- sapply(he_all, function(x){c(0,as.vector(scale(x[2:(length(x)-1)])),0)})
## check for peaks (z-score threshold)
sapply(zscores_he_all, function(x){sum(x > 5)})


#### check uORFs
# on subset: get utr5 lengths, subset
gr_exons <- gr_exons[names(gr_exons) %in% names(gr_cds)]
utr5len <- txLengths[names(gr_cds),]
utr5len <- utr5len[order(rownames(utr5len)),]$utr5_len

leaders <- mapply(function(x,y){x[1:y]}, gr_exons, utr5len)
sapply(leaders, function(x){max(x$ribo_cov)})

#zscore_leaders <- sapply(leaders, function(x){as.vector(scale(x$ribo_cov))})

gr_exons_all <- gr_exons_all[names(gr_exons_all) %in% names(gr_cds_all)]
utr5len_all <- txLengths[names(gr_cds_all),]
utr5len_all <- utr5len_all[order(rownames(utr5len_all)),]$utr5_len

leaders_all <- mapply(function(x,y){x[1:y]}, gr_exons_all, utr5len_all)
sapply(leaders_all, function(x){max(x$ribo_cov)})

###############################
## compare degraded to random non-degraded control
degraded_cds <- gr_cds 
degraded_utr5 <- leaders

non_degraded <- names(gr_cds_all)[!(names(gr_cds_all) %in% names(gr_cds))]
nondeg <- sample(non_degraded, length(degraded_cds))

nondeg_cds <- gr_cds_all[names(gr_cds_all) %in% nondeg]
nondeg_utr5 <- leaders_all[names(leaders_all) %in% nondeg]

deg_avgcov <- sapply(degraded_cds, function(x){sum(x$ribo_cov)/length(x)})
nondeg_avgcov <- sapply(nondeg_cds, function(x){sum(x$ribo_cov)/length(x)})
d <- data.frame(deg = deg_avgcov, nondeg = nondeg_avgcov)
d <- melt(d)
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + scale_x_log10()

zs_deg_cds <- sapply(degraded_cds, function(x){c(0,0,0,as.vector(scale(x$ribo_cov[4:(length(x)-3)])),0,0,0)})
zs_nondeg_cds <- sapply(nondeg_cds, function(x){c(0,0,0,as.vector(scale(x$ribo_cov[4:(length(x)-3)])),0,0,0)})

sum(sapply(zs_deg_cds, function(x){sum(x > 8)}) > 0)
sum(sapply(zs_nondeg_cds, function(x){sum(x > 8)}) > 0)

##### he
deg_he <- zs_deg_cds[names(zs_deg_cds) %in% names(he_all)]
nondeg_he <- zs_nondeg_cds[names(zs_nondeg_cds) %in% names(he_all)]
nondeg_he <- nondeg_he[sample(names(nondeg_he), length(deg_he))]

sum(sapply(deg_he, function(x){sum(x > 8)}) > 0)
sum(sapply(nondeg_he, function(x){sum(x > 8)}) > 0)

#### check patterns around peaks (logo, fld)

