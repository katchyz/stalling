### in silico structure
library(seqinr)
library(GenomicFeatures)
library(data.table)
library(ggplot2)

# get all genes with peaks (longest tx)
# subset fasta: those tx, -30 from start, +30 from end
# save fasta

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
gtf_path <- "../../DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("../../DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "../../DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "../../DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "../../DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "../../DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))

structure_path <- "../../DATA/structure"

###########
### css ###
###########
############ conserved, human
org <- "human"
css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
consensus <- css[!is.na(css$h_ss),]
consensus$seqnames <- consensus$h_tx
consensus$ss <- consensus$h_ss

fasta_cdna <- read.fasta("../../DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
fasta_cdna <- fasta_cdna[names(fasta_cdna) %in% unique(consensus$seqnames)]
# write.fasta(sequences = fasta_cdna, names = names(fasta_cdna), file.out = c(file.path(structure_path, paste0(org, "_conserved.fa"))))

txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]
consensus$cds_len <- sapply(consensus$seqnames, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

mfe <- read.csv(file = c(file.path(structure_path, paste0(org, "_conserved_mfe.txt"))), sep = "\t", header = F)
setDT(mfe)
colnames(mfe) <- c("tx", "mfe")
mfe <- split(mfe, mfe$tx)

##### subset not explained by sequence
load(file = "../../DATA/structure/rest_not_explained_by_AA.Rsave")
rest <- sapply(REST, function(x){strsplit(x, split = "_")[[1]][1]})
mfe <- mfe[names(mfe) %in% rest] # 602 tx

setDT(consensus)
consensus <- consensus[consensus$ss > 30]
consensus <- consensus[consensus$seqnames %in% names(mfe)]
consensus <- consensus[consensus$ss+50 < consensus$cds_len-45]
print(org)
print(nrow(consensus)) # 907

##
consensus$tx_aa <- paste0(consensus$h_tx, "_", (floor(consensus$h_ss/3)+1))
consensus <- consensus[(consensus$tx_aa %in% REST)] # 627

meta <- rep(0,91)
for (j in 1:nrow(consensus)) {
  tx <- as.character(consensus[j,]$seqnames)
  ss <- consensus[j,]$ss
  e <- mfe[[tx]]$mfe[(ss+5-30):(ss+5+60)]
  meta <- meta + e
}

dt <- data.table(sum_mfe = meta, mean_mfe = meta / nrow(consensus), ss_dist = c(-30:60))
dt$position <- rep("stall", nrow(dt))
#ggplot(dt, aes(x = ss_dist, y = mean_mfe)) + geom_line() + ylim(-12,-6) + theme_classic()
#ggsave(c(file.path(structure_path, paste0(org, ".png"))))

### get random position
utx <- unique(as.character(consensus$seqnames))
random_pool <- list()
for (i in 1:length(utx)) {
  tx <- utx[i]
  cons <- consensus[consensus$seqnames == tx]
  cds_len <- cons[1,]$cds_len
  around_ss <- c()
  for (ss in cons$ss) {
    around_ss <- c(around_ss, c((ss-30):(ss+60)))
  }
  if (org == "yeast") {
    pool <- c(30:(cds_len-120))[!c(30:(cds_len-120)) %in% around_ss]
  } else {
    pool <- c(30:(cds_len-120))[!c(30:(cds_len-120)) %in% around_ss]
  }
  pool <- pool[pool %% 3 == 0]
  random_pool[[tx]] <- pool
}

random_pool <- random_pool[sapply(random_pool, function(x){length(x)}) > 0]
#consensus <- consensus[!consensus$seqnames == names(random_pool[sapply(random_pool, function(x){length(x)}) == 0])]
consensus <- consensus[consensus$seqnames %in% names(random_pool)]

###################
## multiple random sampling

random1000 <- data.table(ss_dist = c(-30:60))

for (i in 1:10000) {
  consensus$random1 <- sapply(as.character(consensus$seqnames), function(x){sample(random_pool[[x]],1)})
  meta <- rep(0,91)
  for (j in 1:nrow(consensus)) {
    tx <- as.character(consensus[j,]$seqnames)
    random <- consensus[j,]$random1
    e <- mfe[[tx]]$mfe[(random+5-30):(random+5+60)]
    meta <- meta + e
  }
  mean_meta <- meta / nrow(consensus)
  
  random1000 <- cbind(random1000, mean_meta)

  print(i)
}

# save(random1000, file = "../../DATA/structure/random1000.Rsave")

# make permuted aa control
# load(file = "../../DATA/structure/random1000.Rsave")


dt_r2 <- data.table(sum_mfe = rep(NA, 91), mean_mfe = apply(random1000[,2:10001], 1, mean), ss_dist = c(-30:60))
dt_r2$position <- rep("control", nrow(dt_r2))

dt <- rbind(dt, dt_r2)

ggplot(dt, aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("black", "#615388")) +
  scale_linetype_manual(values = c("dashed", "solid")) + theme_classic()

ggsave(c(file.path(structure_path, "FIXED_CSSs.png")))

dt_permuted$position <- 'permuted'
dt <- rbind(dt, dt_permuted)
ggplot(dt, aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("black", "gray","#615388")) +
  scale_linetype_manual(values = c("dashed", "dashed","solid")) + theme_classic()

dt_p$position <- 'permuted'
dt_final <- rbind(dt[1:182,], dt_permuted)
ggplot(dt_final, aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("black", "gray","#615388")) +
  scale_linetype_manual(values = c("dashed", "dashed","solid")) + theme_classic()

ggsave(c(file.path(structure_path, "CSSs_with_permuted.png")))

consensus$random1 <- sapply(as.character(consensus$seqnames), function(x){sample(random_pool[[x]],1)})
consensus$random2 <- sapply(as.character(consensus$seqnames), function(x){sample(random_pool[[x]],1)})

print(org)
print(nrow(consensus))

meta <- rep(0,91)
for (j in 1:nrow(consensus)) {
  tx <- as.character(consensus[j,]$seqnames)
  random <- consensus[j,]$random1
  e <- mfe[[tx]]$mfe[(random+5-30):(random+5+60)]
  meta <- meta + e
}

dt_r1 <- data.table(sum_mfe = meta, mean_mfe = meta / nrow(consensus), ss_dist = c(-30:60))
dt_r1$position <- rep("random1", nrow(dt_r1))

dt <- rbind(dt, dt_r1)

meta <- rep(0,91)
for (j in 1:nrow(consensus)) {
  tx <- as.character(consensus[j,]$seqnames)
  random <- consensus[j,]$random2
  e <- mfe[[tx]]$mfe[(random+5-30):(random+5+60)]
  meta <- meta + e
}

dt_r2 <- data.table(sum_mfe = meta, mean_mfe = meta / nrow(consensus), ss_dist = c(-30:60))
dt_r2$position <- rep("random2", nrow(dt_r2))

dt <- rbind(dt, dt_r2)

ggplot(dt, aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()

ggsave(c(file.path(structure_path, "CSSs.png")))


for (j in 1:nrow(consensus)) {
  tx <- as.character(consensus[j,]$seqnames)
  ss <- consensus[j,]$ss
  print(fasta_fold[[tx]][(ss+5):(ss+5+2)])
}

dt_q5 <- data.table(sum_mfe = rep(NA, 91), mean_mfe = apply(random1000[,2:10001], 1, function(x) quantile(x, 0.05)), ss_dist = c(-30:60))
dt_q5$position <- rep("q5", nrow(dt_q5))
dt_q95 <- data.table(sum_mfe = rep(NA, 91), mean_mfe = apply(random1000[,2:10001], 1, function(x) quantile(x, 0.95)), ss_dist = c(-30:60))
dt_q95$position <- rep("q95", nrow(dt_q95))

dt_quan <- dt_final
dt_quan <- rbind(dt_quan, dt_q5)
dt_quan <- rbind(dt_quan, dt_q95)

dt_quan$qmin <- rep(dt_quan[dt_quan$position == "q5"]$mean_mfe, 5)
dt_quan$qmax <- rep(dt_quan[dt_quan$position == "q95"]$mean_mfe, 5)

ggplot(dt_quan[1:273], aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("black", "gray", "#615388")) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid")) + theme_classic() +
  geom_ribbon(data=subset(dt_quan, position == "stall"), aes(ymin=qmin, ymax=qmax),fill="gray", alpha=0.2)

ggsave(c(file.path(structure_path, "CSSs_with_quan.png")))
save(dt_quan, file = "../../DATA/structure/dt_quan.Rsave")

dt_notp <- dt_quan[dt_quan$position != "permuted"]
dt_notp <- rbind(dt_notp[,1:4], dt_p)
dt_all <- rbind(dt_notp, dt_q5_perm) 
dt_all <- rbind(dt_all, dt_q95_perm) 
dt_all$qmin <- rep(dt_all[dt_all$position == "q5"]$mean_mfe, 7)
dt_all$qmax <- rep(dt_all[dt_all$position == "q95"]$mean_mfe, 7)
dt_all$qminp <- rep(dt_all[dt_all$position == "q5_perm"]$mean_mfe, 7)
dt_all$qmaxp <- rep(dt_all[dt_all$position == "q95_perm"]$mean_mfe, 7)

ggplot(dt_all[dt_all$position %in% c("stall", "control", "permuted")], aes(x = ss_dist, y = mean_mfe)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("black", "gray", "#615388")) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid")) + theme_classic() +
  geom_ribbon(data=subset(dt_all, position == "stall"), aes(ymin=qmin, ymax=qmax),fill="black", alpha=0.2) +
  geom_ribbon(data=subset(dt_all, position == "stall"), aes(ymin=qminp, ymax=qmaxp),fill="gray", alpha=0.2)

ggsave(c(file.path(structure_path, "CSSs_with_distributions.png")))

