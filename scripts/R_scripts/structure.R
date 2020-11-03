### structure
library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(seqinr)

# RNAfold-ed
fasta_vienna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures_only.txt")
fasta_vienna <- fasta_vienna[names(fasta_vienna) %in% sam]
fasta_vienna <- sapply(fasta_vienna, function(x){recode(getSequence(x),"'.'=1;'('=0;')'=0")})

# DMS-seq
dmsseq_path <- "/Volumes/USELESS/DATA/DMSseq"

dms_2h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_2h_2 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_b1_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_2h_3 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_b2_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")

dms_4h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_sphere_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_4h_2 <- read.csv(file = c(file.path(dmsseq_path, "WT_sphere_b1_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_4h_3 <- read.csv(file = c(file.path(dmsseq_path, "WT_sphere_b2_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")

dms_6h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_shield_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_6h_2 <- read.csv(file = c(file.path(dmsseq_path, "WT_shield_b1_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
dms_6h_3 <- read.csv(file = c(file.path(dmsseq_path, "WT_shield_b2_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")

# convert to granges (for faster loading and analysis)
make_granges <- function(x) {
  GRanges(seqnames=x[1],
          IRanges(start=1:x[4], width=1),
          accessibility=as.numeric(strsplit(as.character(x[10]), split = ",")[[1]]),
          seq = strsplit(as.character(x[11]), split = ""))
}

# 2h
gr_2h_1 <- GRangesList(apply(dms_2h_1, 1, make_granges))
names(gr_2h_1) <- sapply(gr_2h_1, function(x){as.character(runValue(seqnames(x)))})
save(gr_2h_1, file = c(file.path(dmsseq_path, "gr_2h_1.Rsave")))

gr_2h_2 <- GRangesList(apply(dms_2h_2, 1, make_granges))
names(gr_2h_2) <- sapply(gr_2h_2, function(x){as.character(runValue(seqnames(x)))})
save(gr_2h_2, file = c(file.path(dmsseq_path, "gr_2h_2.Rsave")))

gr_2h_3 <- GRangesList(apply(dms_2h_3, 1, make_granges))
names(gr_2h_3) <- sapply(gr_2h_3, function(x){as.character(runValue(seqnames(x)))})
save(gr_2h_3, file = c(file.path(dmsseq_path, "gr_2h_3.Rsave")))

# 4h
gr_4h_1 <- GRangesList(apply(dms_4h_1, 1, make_granges))
names(gr_4h_1) <- sapply(gr_4h_1, function(x){as.character(runValue(seqnames(x)))})
save(gr_4h_1, file = c(file.path(dmsseq_path, "gr_4h_1.Rsave")))

gr_4h_2 <- GRangesList(apply(dms_4h_2, 1, make_granges))
names(gr_4h_2) <- sapply(gr_4h_2, function(x){as.character(runValue(seqnames(x)))})
save(gr_4h_2, file = c(file.path(dmsseq_path, "gr_4h_2.Rsave")))

gr_4h_3 <- GRangesList(apply(dms_4h_3, 1, make_granges))
names(gr_4h_3) <- sapply(gr_4h_3, function(x){as.character(runValue(seqnames(x)))})
save(gr_4h_3, file = c(file.path(dmsseq_path, "gr_4h_3.Rsave")))

# 6h
gr_6h_1 <- GRangesList(apply(dms_6h_1, 1, make_granges))
names(gr_6h_1) <- sapply(gr_6h_1, function(x){as.character(runValue(seqnames(x)))})
save(gr_6h_1, file = c(file.path(dmsseq_path, "gr_6h_1.Rsave")))

gr_6h_2 <- GRangesList(apply(dms_6h_2, 1, make_granges))
names(gr_6h_2) <- sapply(gr_6h_2, function(x){as.character(runValue(seqnames(x)))})
save(gr_6h_2, file = c(file.path(dmsseq_path, "gr_6h_2.Rsave")))

gr_6h_3 <- GRangesList(apply(dms_6h_3, 1, make_granges))
names(gr_6h_3) <- sapply(gr_6h_3, function(x){as.character(runValue(seqnames(x)))})
save(gr_6h_3, file = c(file.path(dmsseq_path, "gr_6h_3.Rsave")))

#######################################################################################
### metaplots

# load DMS-seq data
dmsseq_path <- "/Volumes/USELESS/DATA/DMSseq"
libs <- list.files(path = dmsseq_path, pattern = "\\.Rsave$")
for (l in 1:length(libs)) {
  load(c(file.path(dmsseq_path, libs[l])))
}

# get zebrafish stall sites
zebrafihs_path <- "/Volumes/USELESS/STALLING/zebrafish_regulatory_ss"
load(file = c(file.path(zebrafihs_path, "consensus_2h.Rsave")))
load(file = c(file.path(zebrafihs_path, "consensus_4h.Rsave")))
load(file = c(file.path(zebrafihs_path, "consensus_6h.Rsave")))

# transcript db
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.81_chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
## longest tx
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name

### use only 1st replicate (most genes in there)
sum(consensus_2h$seqnames %in% names(gr_2h_1)) # 73
sum(consensus_4h$seqnames %in% names(gr_4h_1)) # 49
sum(consensus_6h$seqnames %in% names(gr_6h_1)) # 51

# get position on CDS
consensus_2h <- consensus_2h[consensus_2h$seqnames %in% names(gr_2h_1)]
consensus_2h$utr5_len <- sapply(consensus_2h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
consensus_2h$ss <- consensus_2h$start - consensus_2h$utr5_len
consensus_2h <- consensus_2h[consensus_2h$ss > 15,] # 66
consensus_2h$cds_len <- sapply(consensus_2h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$cds_len})
consensus_2h$gene_id <- sapply(consensus_2h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$gene_id})

consensus_4h <- consensus_4h[consensus_4h$seqnames %in% names(gr_4h_1)]
consensus_4h$utr5_len <- sapply(consensus_4h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
consensus_4h$ss <- consensus_4h$start - consensus_4h$utr5_len
consensus_4h <- consensus_4h[consensus_4h$ss > 15,] # 43
consensus_4h$cds_len <- sapply(consensus_4h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$cds_len})
consensus_4h$gene_id <- sapply(consensus_4h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$gene_id})

consensus_6h <- consensus_6h[consensus_6h$seqnames %in% names(gr_6h_1)]
consensus_6h$utr5_len <- sapply(consensus_6h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
consensus_6h$ss <- consensus_6h$start - consensus_6h$utr5_len
consensus_6h <- consensus_6h[consensus_6h$ss > 15,] # 40
consensus_6h$cds_len <- sapply(consensus_6h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$cds_len})
consensus_6h$gene_id <- sapply(consensus_6h$seqnames, function(x){txLengths[txLengths$tx_name == x,]$gene_id})

### different transcripts/genes at different stages!

### metaplots (use "start", as DMS-seq is in CDNA coordinates)

# 2h
meta_2h <- rep(0,115)
ac_2h <- rep(0,115)
for (i in 1:nrow(consensus_2h)) {
  tx <- as.character(consensus_2h[i,]$seqnames)
  ss <- consensus_2h[i,]$start
  
  meta_ss <- gr_2h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ac_2h <- ac_2h + as.integer(!is.na(meta_ss))
  meta_ss <- replace(meta_ss, is.na(meta_ss), 0)
  
  meta_2h <- meta_2h + meta_ss
}

scale <- c(-15:99)
meta_2h <- data.table(accessibility = meta_2h, dist_ss = scale, ac = ac_2h / nrow(consensus_2h))
p2h <- ggplot(meta, aes(x = dist_ss, y = accessibility*ac)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/dms_2h.png", plot = p2h)

# 4h
meta_4h <- rep(0,115)
ac_4h <- rep(0,115)
for (i in 1:nrow(consensus_4h)) {
  tx <- as.character(consensus_4h[i,]$seqnames)
  ss <- consensus_4h[i,]$start
  
  meta_ss <- gr_4h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ac_4h <- ac_4h + as.integer(!is.na(meta_ss))
  meta_ss <- replace(meta_ss, is.na(meta_ss), 0)
  
  meta_4h <- meta_4h + meta_ss
}

scale <- c(-15:99)
meta_4h <- data.table(accessibility = meta_4h, dist_ss = scale, ac = ac_4h / nrow(consensus_4h))
p4h <- ggplot(meta_4h, aes(x = dist_ss, y = accessibility)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/dms_4h.png", plot = p4h)

# 6h
meta_6h <- rep(0,115)
ac_6h <- rep(0,115)
for (i in 1:nrow(consensus_6h)) {
  tx <- as.character(consensus_6h[i,]$seqnames)
  ss <- consensus_6h[i,]$start
  
  meta_ss <- gr_6h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ac_6h <- ac_6h + as.integer(!is.na(meta_ss))
  meta_ss <- replace(meta_ss, is.na(meta_ss), 0)
  
  meta_6h <- meta_6h + meta_ss
}

scale <- c(-15:99)
meta_6h <- data.table(seqnames = consensus_6h$seqnames, accessibility = meta_6h, dist_ss = scale, ac = ac_6h / nrow(consensus_6h))
p6h <- ggplot(meta, aes(x = dist_ss, y = accessibility)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/dms_6h.png", plot = p6h)

#### smoothing
meta_2h$smoothed <- smooth(meta_2h$accessibility*meta_2h$ac)
meta_4h$smoothed <- smooth(meta_4h$accessibility*meta_4h$ac)
meta_6h$smoothed <- smooth(meta_6h$accessibility*meta_6h$ac)

s2h <- ggplot(meta_2h, aes(x = dist_ss, y = smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/smoothed_2h.png", plot = s2h)

s4h <- ggplot(meta_4h, aes(x = dist_ss, y = smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/smoothed_4h.png", plot = s4h)

s6h <- ggplot(meta_6h, aes(x = dist_ss, y = smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/smoothed_6h.png", plot = s6h)

####################################################################################################
##### get CONTROL (DMS, untreated)
dmsseq_path <- "/Volumes/USELESS/DATA/DMSseq"

ctr_2h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_DMS_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
ctr_2h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_DMS_B1_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")
ctr_2h_1 <- read.csv(file = c(file.path(dmsseq_path, "WT_64c_DMS_B2_accessibility_profiles_78_protein_coding_genes.csv")),
                     header = TRUE, sep = "\t")


ctr_gr_2h_1 <- GRangesList(apply(ctr_2h_1, 1, make_granges))
names(ctr_gr_2h_1) <- sapply(ctr_gr_2h_1, function(x){as.character(runValue(seqnames(x)))})
save(ctr_gr_2h_1, file = c(file.path(dmsseq_path, "ctr_gr_2h_1.Rsave")))


# control 2h
ctr_meta_2h <- rep(0,115)
ctr_ac_2h <- rep(0,115)
for (i in 1:nrow(consensus_2h)) {
  tx <- as.character(consensus_2h[i,]$seqnames)
  ss <- consensus_2h[i,]$start
  
  ctr_meta_ss <- ctr_gr_2h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ctr_ac_2h <- ctr_ac_2h + as.integer(!is.na(ctr_meta_ss))
  ctr_meta_ss <- replace(ctr_meta_ss, is.na(ctr_meta_ss), 0)
  
  ctr_meta_2h <- ctr_meta_2h + ctr_meta_ss
}

scale <- c(-15:99)
ctr_meta_2h <- data.table(accessibility = ctr_meta_2h, dist_ss = scale, ac = ctr_ac_2h / nrow(consensus_2h))
c2h <- ggplot(ctr_meta_2h, aes(x = dist_ss, y = accessibility*ac)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/ctr_2h.png", plot = c2h)

ctr_meta_2h$smoothed <- smooth(ctr_meta_2h$accessibility*ctr_meta_2h$ac)
sc2h <- ggplot(ctr_meta_2h, aes(x = dist_ss, y = smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/ctr_smoothed_2h.png", plot = sc2h)

############## plot DMS-treated (2/4/6h) minus DMS-untreated (2h)

meta_2h$ctr <- ctr_meta_2h$accessibility
meta_2h$ctr_smoothed <- ctr_meta_2h$smoothed
ggplot(meta_2h, aes(x = dist_ss, y = (accessibility*ac)-(ctr*ac))) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_2h.png")
ggplot(meta_2h, aes(x = dist_ss, y = smoothed-ctr_smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_smoothed_2h.png")

# control 4h
consensus_4h <- consensus_4h[consensus_4h$seqnames %in% names(ctr_gr_2h_1),]

# 4h
meta_4h <- rep(0,115)
ac_4h <- rep(0,115)
for (i in 1:nrow(consensus_4h)) {
  tx <- as.character(consensus_4h[i,]$seqnames)
  ss <- consensus_4h[i,]$start
  
  meta_ss <- gr_4h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ac_4h <- ac_4h + as.integer(!is.na(meta_ss))
  meta_ss <- replace(meta_ss, is.na(meta_ss), 0)
  
  meta_4h <- meta_4h + meta_ss
}

scale <- c(-15:99)
meta_4h <- data.table(accessibility = meta_4h, dist_ss = scale, ac = ac_4h / nrow(consensus_4h))
meta_4h$smoothed <- smooth(meta_4h$accessibility*meta_4h$ac)

ctr_meta_4h <- rep(0,115)
ctr_ac_4h <- rep(0,115)
for (i in 1:nrow(consensus_4h)) {
  tx <- as.character(consensus_4h[i,]$seqnames)
  ss <- consensus_4h[i,]$start
  
  ctr_meta_ss <- ctr_gr_2h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ctr_ac_4h <- ctr_ac_4h + as.integer(!is.na(ctr_meta_ss))
  ctr_meta_ss <- replace(ctr_meta_ss, is.na(ctr_meta_ss), 0)
  
  ctr_meta_4h <- ctr_meta_4h + ctr_meta_ss
}
scale <- c(-15:99)
ctr_meta_4h <- data.table(accessibility = ctr_meta_4h, dist_ss = scale, ac = ctr_ac_4h / nrow(consensus_4h))
ctr_meta_4h$smoothed <- smooth(ctr_meta_4h$accessibility*ctr_meta_4h$ac)

meta_4h$ctr <- ctr_meta_4h$accessibility
meta_4h$ctr_smoothed <- ctr_meta_4h$smoothed

ggplot(meta_4h, aes(x = dist_ss, y = (accessibility*ac)-(ctr*ac))) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_4h.png")
ggplot(meta_4h, aes(x = dist_ss, y = smoothed-ctr_smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_smoothed_4h.png")

#############
# control 6h
consensus_6h <- consensus_6h[consensus_6h$seqnames %in% names(ctr_gr_2h_1),]

# 6h
meta_6h <- rep(0,115)
ac_6h <- rep(0,115)
for (i in 1:nrow(consensus_6h)) {
  tx <- as.character(consensus_6h[i,]$seqnames)
  ss <- consensus_6h[i,]$start
  
  meta_ss <- gr_6h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ac_6h <- ac_6h + as.integer(!is.na(meta_ss))
  meta_ss <- replace(meta_ss, is.na(meta_ss), 0)
  
  meta_6h <- meta_6h + meta_ss
}

scale <- c(-15:99)
meta_6h <- data.table(accessibility = meta_6h, dist_ss = scale, ac = ac_6h / nrow(consensus_6h))
meta_6h$smoothed <- smooth(meta_6h$accessibility*meta_6h$ac)

ctr_meta_6h <- rep(0,115)
ctr_ac_6h <- rep(0,115)
for (i in 1:nrow(consensus_6h)) {
  tx <- as.character(consensus_6h[i,]$seqnames)
  ss <- consensus_6h[i,]$start
  
  ctr_meta_ss <- ctr_gr_2h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  ctr_ac_6h <- ctr_ac_6h + as.integer(!is.na(ctr_meta_ss))
  ctr_meta_ss <- replace(ctr_meta_ss, is.na(ctr_meta_ss), 0)
  
  ctr_meta_6h <- ctr_meta_6h + ctr_meta_ss
}
scale <- c(-15:99)
ctr_meta_6h <- data.table(accessibility = ctr_meta_6h, dist_ss = scale, ac = ctr_ac_6h / nrow(consensus_6h))
ctr_meta_6h$smoothed <- smooth(ctr_meta_6h$accessibility*ctr_meta_6h$ac)

meta_6h$ctr <- ctr_meta_6h$accessibility
meta_6h$ctr_smoothed <- ctr_meta_6h$smoothed

ggplot(meta_6h, aes(x = dist_ss, y = (accessibility*ac)-(ctr*ac))) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_6h.png")
ggplot(meta_6h, aes(x = dist_ss, y = smoothed-ctr_smoothed)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/STALLING/PLOTS/structure/diff_smoothed_6h.png")


##### normalize
dupa_2h <- rep(0,115)
dac_2h <- rep(0,115)
for (i in 1:nrow(consensus_2h)) {
  tx <- as.character(consensus_2h[i,]$seqnames)
  ss <- consensus_2h[i,]$start
  
  dupa_ss <- gr_2h_1[[tx]]$accessibility[(ss-15):(ss+99)]
  dac_2h <- dac_2h + as.integer(!is.na(dupa_ss))
  dupa_ss <- replace(dupa_ss, is.na(dupa_ss), 0)
  dupa_ss <- dupa_ss / sum(dupa_ss)
  dupa_2h <- dupa_2h + dupa_ss
}

scale <- c(-15:99)
dupa_2h <- data.table(accessibility = dupa_2h, dist_ss = scale, ac = dac_2h / nrow(consensus_2h))
ggplot(dupa_2h, aes(x = dist_ss, y = accessibility*ac)) + geom_bar(stat = "identity")
dupa_2h$smoothed <- smooth(dupa_2h$accessibility*dupa_2h$ac)
ggplot(dupa_2h, aes(x = dist_ss, y = smoothed)) + geom_bar(stat = "identity")


################################################

### cldnd "ENSDART00000021620" (is in consensus_2h) (is in gr_2h and gr_4h)
cldnd <- data.table(acc_2h_1 = gr_2h_1[["ENSDART00000021620"]]$accessibility,
                    acc_2h_2 = gr_2h_2[["ENSDART00000021620"]]$accessibility,
                    acc_2h_3 = gr_2h_3[["ENSDART00000021620"]]$accessibility,
                    acc_4h_1 = gr_4h_1[["ENSDART00000021620"]]$accessibility,
                    acc_4h_2 = gr_4h_2[["ENSDART00000021620"]]$accessibility,
                    acc_4h_3 = gr_4h_3[["ENSDART00000021620"]]$accessibility,
                    scale = c(1:777))
cldnd[is.na(cldnd)] <- 0

ggplot(cldnd[100:175], aes(x = scale, y = acc_2h_1)) + geom_bar(stat = "identity")
ggplot(cldnd[100:175], aes(x = scale, y = acc_4h_1)) + geom_bar(stat = "identity")
ggplot(cldnd[100:175], aes(x = scale, y = acc_4h_1-acc_2h_1)) + geom_bar(stat = "identity")

### SHAPES
SHAPE_cell256 <- get(load(file = "/Volumes/USELESS/DATA/SHAPES/SE/norm2_cell256_invivo.Rsave"))
SHAPE_oblong <- get(load(file = "/Volumes/USELESS/DATA/SHAPES/SE/norm2_oblong_invivo.Rsave"))

cldnd_SHAPE <- data.table(cell256_treated = SHAPE_cell256[["ENSDART00000021620"]]$TC.treated,
                          cell256_control = SHAPE_cell256[["ENSDART00000021620"]]$TC.control,
                          cell256_log2ratio = SHAPE_cell256[["ENSDART00000021620"]]$log2ratio,
                          oblong_treated = SHAPE_oblong[["ENSDART00000021620"]]$TC.treated,
                          oblong_control = SHAPE_oblong[["ENSDART00000021620"]]$TC.control,
                          oblong_log2ratio = SHAPE_oblong[["ENSDART00000021620"]]$log2ratio,
                          scale = c(1:777))

ggplot(cldnd_SHAPE[100:175], aes(x = scale, y = cell256_log2ratio-oblong_log2ratio)) + geom_bar(stat = "identity")

