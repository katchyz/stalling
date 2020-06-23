### zebrafish - any different stall sites?
# (1) get peaks at each stage (in min.2/both libraries)
# (2) check level of expression of these genes (tables / ribo-seq / rna-seq / mishima)
# (3) see if any peak clearly comes up in one stage, but not others

##################################################################
#### (1) ####
library(data.table)
library(dplyr)
library(plyr)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
zebrafihs_path <- "/Volumes/USELESS/STALLING/zebrafish_regulatory_ss"

org <- "zebrafish"

libs_2h <- c("zebrafish_2h.Rsave", "zebrafish_beaudoin_2h.Rsave", "zebrafish_subtelny_2h.Rsave")
libs_4h <- c("zebrafish_dome.Rsave", "zebrafish_subtelny_4h.Rsave")
libs_6h <- c("zebrafish_shield.Rsave", "zebrafish_subtelny_6h.Rsave")

### longest tx
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# strand
strand <- transcripts(txdb)
strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
## longest tx
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name

###
dt_2h <- data.table(seqnames=c("transcript"), start=c(0))
for (j in 1:length(libs_2h)){
  load(file = c(file.path(peaks_path, libs_2h[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt_2h <- rbind(dt_2h, peaks)
}
consensus_2h <- plyr::count(dt_2h, vars = c("seqnames", "start"))
consensus_2h <- consensus_2h[consensus_2h$freq > 1,]
setDT(consensus_2h)
consensus_2h <- consensus_2h[,1:2]
save(consensus_2h, file = c(file.path(zebrafihs_path, "consensus_2h.Rsave")))
###
dt_4h <- data.table(seqnames=c("transcript"), start=c(0))
for (j in 1:length(libs_4h)){
  load(file = c(file.path(peaks_path, libs_4h[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt_4h <- rbind(dt_4h, peaks)
}
consensus_4h <- plyr::count(dt_4h, vars = c("seqnames", "start"))
consensus_4h <- consensus_4h[consensus_4h$freq > 1,]
setDT(consensus_4h)
consensus_4h <- consensus_4h[,1:2]
save(consensus_4h, file = c(file.path(zebrafihs_path, "consensus_4h.Rsave")))
###
dt_6h <- data.table(seqnames=c("transcript"), start=c(0))
for (j in 1:length(libs_6h)){
  load(file = c(file.path(peaks_path, libs_6h[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt_6h <- rbind(dt_6h, peaks)
}
consensus_6h <- plyr::count(dt_6h, vars = c("seqnames", "start"))
consensus_6h <- consensus_6h[consensus_6h$freq > 1,]
setDT(consensus_6h)
consensus_6h <- consensus_6h[,1:2]
save(consensus_6h, file = c(file.path(zebrafihs_path, "consensus_6h.Rsave")))

################################
#### (2,3) ####

mishima <- read.csv(file = "/Volumes/USELESS/STALLING/zebrafish_regulatory_ss/mishima/stable.csv", header = T)
mishima <- mishima[,1:2]
setDT(mishima)

# get gene names of my transcripts
consensus_2h$gene_id <- sapply(consensus_2h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 117
consensus_4h$gene_id <- sapply(consensus_4h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 73
consensus_6h$gene_id <- sapply(consensus_6h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 73

m2h <- consensus_2h[consensus_2h$gene_id %in% mishima$gene_id] # 20
m4h <- consensus_4h[consensus_4h$gene_id %in% mishima$gene_id] # 13
m6h <- consensus_6h[consensus_6h$gene_id %in% mishima$gene_id] # 10

# check present/missing stall sites
dt <- data.table(seqnames=c("transcript"), start=c(0), gene_id=c("gene"))
dt <- rbind(dt, m2h)
dt <- rbind(dt, m4h)
dt <- rbind(dt, m6h)
consensus <- plyr::count(dt, vars = c("seqnames", "start"))


# present in 2: ENSDART00000113351 at 143, ENSDART00000101109 at 186, ENSDART00000114762 at 553, ENSDART00000158367 at 2833/2836
# present in 1: all the rest

# add gene name
m2h$gene <- sapply(m2h$gene_id, function(x){mishima[mishima$gene_id == x,]$gene})
m4h$gene <- sapply(m4h$gene_id, function(x){mishima[mishima$gene_id == x,]$gene})
m6h$gene <- sapply(m6h$gene_id, function(x){mishima[mishima$gene_id == x,]$gene})

save(m2h, file = c(file.path(zebrafihs_path, "stable_ss_2h.Rsave")))
save(m4h, file = c(file.path(zebrafihs_path, "stable_ss_4h.Rsave")))
save(m6h, file = c(file.path(zebrafihs_path, "stable_ss_6h.Rsave")))

# add cds coordinates (amino acid number)
m2h$utr5_len <- sapply(m2h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
m4h$utr5_len <- sapply(m4h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
m6h$utr5_len <- sapply(m6h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})

m2h$ss <- m2h$start - m2h$utr5_len
m4h$ss <- m4h$start - m4h$utr5_len
m6h$ss <- m6h$start - m6h$utr5_len

m2h$aa <- floor(m2h$ss / 3) + 1
m4h$aa <- floor(m4h$ss / 3) + 1
m6h$aa <- floor(m6h$ss / 3) + 1

m2h <- m2h[m2h$ss > 15,]
m4h <- m4h[m4h$ss > 15,]
m6h <- m6h[m6h$ss > 15,]

### make sure the stall sites do not come up in (single) libraries from other stages
# for m2h - peaks 4 & 6h, etc.
u2h <- unique(dt_2h[2:nrow(dt_2h),])
u4h <- unique(dt_4h[2:nrow(dt_4h),])
u6h <- unique(dt_6h[2:nrow(dt_6h),])

m2h_u4h <- plyr::count(rbind(m2h[,1:2], u4h), vars = c("seqnames", "start"))
m2h_u4h[m2h_u4h$freq == 2,]
m2h_u6h <- plyr::count(rbind(m2h[,1:2], u6h), vars = c("seqnames", "start"))
m2h_u6h[m2h_u6h$freq == 2,]

m4h_u2h <- plyr::count(rbind(m4h[,1:2], u2h), vars = c("seqnames", "start"))
m4h_u2h[m4h_u2h$freq == 2,]
m4h_u6h <- plyr::count(rbind(m4h[,1:2], u6h), vars = c("seqnames", "start"))
m4h_u6h[m4h_u6h$freq == 2,]

m6h_u2h <- plyr::count(rbind(m6h[,1:2], u2h), vars = c("seqnames", "start"))
m6h_u2h[m6h_u2h$freq == 2,]
m6h_u4h <- plyr::count(rbind(m6h[,1:2], u4h), vars = c("seqnames", "start"))
m6h_u4h[m6h_u4h$freq == 2,]

######## clearly in 2libs of one stage - and in none of the libraries of at least one another stage
#############################################
# high in 2h, low in 6h
mdecay <- read.csv(file = "/Volumes/USELESS/STALLING/zebrafish_regulatory_ss/mishima/M-decay.csv", header = T)
mdecay <- mdecay[,1:2]
setDT(mdecay)

# high in alpha-amanitin 6h embryos compared to 6h WT (mRNA stabilized in the absence of zygotic translation)
zdecay <- read.csv(file = "/Volumes/USELESS/STALLING/zebrafish_regulatory_ss/mishima/Z-decay.csv", header = T)
zdecay <- zdecay[,1:2]
setDT(zdecay)

load(file = c(file.path(zebrafihs_path, "consensus_2h.Rsave")))
load(file = c(file.path(zebrafihs_path, "consensus_4h.Rsave")))
load(file = c(file.path(zebrafihs_path, "consensus_6h.Rsave")))
# get gene names of my transcripts
consensus_2h$gene_id <- sapply(consensus_2h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 117
consensus_4h$gene_id <- sapply(consensus_4h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 73
consensus_6h$gene_id <- sapply(consensus_6h$seqnames, function(x){txlen[txlen$tx_name == x,]$gene_id}) # 73

m2h <- consensus_2h[consensus_2h$gene_id %in% mdecay$gene_id] # 4
m4h <- consensus_4h[consensus_4h$gene_id %in% mdecay$gene_id] # 2
m6h <- consensus_6h[consensus_6h$gene_id %in% mdecay$gene_id] # 0

z2h <- consensus_2h[consensus_2h$gene_id %in% zdecay$gene_id] # 5
z4h <- consensus_4h[consensus_4h$gene_id %in% zdecay$gene_id] # 3
z6h <- consensus_6h[consensus_6h$gene_id %in% zdecay$gene_id] # 0

# add cds coordinates (amino acid number)
m2h$utr5_len <- sapply(m2h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
m4h$utr5_len <- sapply(m4h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
m6h$utr5_len <- sapply(m6h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})

m2h$ss <- m2h$start - m2h$utr5_len
m4h$ss <- m4h$start - m4h$utr5_len
m6h$ss <- m6h$start - m6h$utr5_len

m2h$aa <- floor(m2h$ss / 3) + 1
m4h$aa <- floor(m4h$ss / 3) + 1
m6h$aa <- floor(m6h$ss / 3) + 1

m2h <- m2h[m2h$ss > 15,]
m4h <- m4h[m4h$ss > 15,]
m6h <- m6h[m6h$ss > 15,]

# add cds coordinates (amino acid number)
z2h$utr5_len <- sapply(z2h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
z4h$utr5_len <- sapply(z4h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})
z6h$utr5_len <- sapply(z6h$seqnames, function(x){txlen[txlen$tx_name == x,]$utr5_len})

z2h$ss <- z2h$start - z2h$utr5_len
z4h$ss <- z4h$start - z4h$utr5_len
z6h$ss <- z6h$start - z6h$utr5_len

z2h$aa <- floor(z2h$ss / 3) + 1
z4h$aa <- floor(z4h$ss / 3) + 1
z6h$aa <- floor(z6h$ss / 3) + 1

z2h <- z2h[z2h$ss > 15,]
z4h <- z4h[z4h$ss > 15,]
z6h <- z6h[z6h$ss > 15,]

### add gene name
m2h$gene <- sapply(m2h$gene_id, function(x){mdecay[mdecay$gene_id == x,]$gene})
m4h$gene <- sapply(m4h$gene_id, function(x){mdecay[mdecay$gene_id == x,]$gene})
m6h$gene <- sapply(m6h$gene_id, function(x){mdecay[mdecay$gene_id == x,]$gene})

z2h$gene <- sapply(z2h$gene_id, function(x){zdecay[zdecay$gene_id == x,]$gene})
z4h$gene <- sapply(z4h$gene_id, function(x){zdecay[zdecay$gene_id == x,]$gene})
z6h$gene <- sapply(z6h$gene_id, function(x){zdecay[zdecay$gene_id == x,]$gene})

