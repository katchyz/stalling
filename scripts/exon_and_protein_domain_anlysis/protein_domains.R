### protein domain boundaries
library(data.table)
library(ggplot2)

cath <- read.table("../../DATA/domains/cath_human.tsv", header = T, sep = "\t")


ss_uniprot <- read.table("./../DATA/domains/uniprot_gene.tab",
                         header = T, sep = "\t")
ss_uniprot$Gene.names <- sapply(as.character(ss_uniprot$Gene.names), function(x){strsplit(x, split = " ")[[1]][1]})

cath <- cath[cath$UNIPROT_ACC %in% ss_uniprot$Entry,]
setDT(cath)

##
dt <- data.table(uniprot = cath$UNIPROT_ACC, boundaries = cath$BOUNDARIES)
dt <- unique(dt)

dt$boundaries <- as.character(dt$boundaries)
multb <- dt[grepl(",", dt$boundaries)]
dt <- dt[!grepl(",", dt$boundaries)]

singb <- data.table(uniprot = c(), boundaries = c())
for (i in 1:nrow(multb)) {
  up <- multb[i,]$uniprot
  boundaries <- multb[i,]$boundaries
  boundaries <- strsplit(boundaries, split = ",")[[1]]
  for (j in 1:length(boundaries)) {
    singb <- rbind(singb, data.table(uniprot = up, boundaries = boundaries[j]))
  }
}

dt <- rbind(dt, singb)
dt <- dt[order(dt$uniprot)]
dt <- unique(dt)
dt$gene <- sapply(as.character(dt$uniprot), function(x){ss_uniprot[ss_uniprot$Entry == x,]$Gene.names})
dt$boundary_start <- as.numeric(sapply(dt$boundaries, function(x){strsplit(x, split = "-")[[1]][1]}))
dt$boundary_end <- as.numeric(sapply(dt$boundaries, function(x){strsplit(x, split = "-")[[1]][2]}))

css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
human_ss <- css[!is.na(css$h_ss),][,10:12]
human_ss$aa <- floor(human_ss$h_ss / 3) + 1

dist_to_bound_end <- data.table(gene = c(), dist = c())
for (i in 1:nrow(human_ss)) {
  geneN <- human_ss[i,]$gene
  aa <- human_ss[i,]$aa
  be <- dt[dt$gene == geneN]$boundary_end
  dist_be <- data.table(gene = rep(geneN, length(be)), dist = (aa - be))
  dist_to_bound_end <- rbind(dist_to_bound_end, dist_be)
  print(i)
}

dist_pos <- dist_to_bound_end[dist_to_bound_end$dist > -1]
ggplot(dist_pos, aes(x = dist_pos$dist)) + geom_histogram()

min_dist_to_bound_end <- data.table(gene = c(), dist = c())
for (i in 1:nrow(human_ss)) {
  geneN <- human_ss[i,]$gene
  aa <- human_ss[i,]$aa
  be <- dt[dt$gene == geneN]$boundary_end
  dist_be <- aa - be
  dist_be <- dist_be[dist_be > -1]
  if (length(dist_be) > 0) {
    dist_be <- min(dist_be)
  } else {
    dist_be <- NA
  }
  dist_be <- data.table(gene = geneN, dist = dist_be)
  min_dist_to_bound_end <- rbind(min_dist_to_bound_end, dist_be)
  print(i)
}

ggplot(min_dist_to_bound_end, aes(x = dist)) + geom_histogram(binwidth = 1) + xlim(0,500)
# ggsave("protein_boundaries.png")


######### control
# add random control positions to human_ss
library(GenomicFeatures)
library(reshape2)

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]

human_ss$cds_len <- sapply(human_ss$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
human_ss$cds_len <- floor(human_ss$cds_len / 3) # in amino acids

### CONTROL
utx <- unique(as.character(human_ss$h_tx))
random_pool <- list()
for (i in 1:length(utx)) {
  txn <- utx[i]
  cds_len_aa <- human_ss[human_ss$h_tx == txn,]$cds_len[1]
  around_ss <- c()
  for (ss in human_ss[as.character(human_ss$h_tx) == txn,]$aa) {
    around_ss <- c(around_ss, c((ss-5):(ss+6)))
  }
  pool <- c(6:(cds_len_aa-5))[!c(6:(cds_len_aa-5)) %in% around_ss]
  random_pool[[txn]] <- pool
}

random_pool <- random_pool[sapply(random_pool, function(x){length(x)}) > 0] # all 1729

###### 1000 random
min_dist_to_bound_end <- data.table(ss = rep(0, nrow(human_ss)))
for (i in 1:nrow(human_ss)) {
  geneN <- human_ss[i,]$gene
  aa <- human_ss[i,]$aa
  be <- dt[dt$gene == geneN]$boundary_end
  # ss
  dist_be <- aa - be
  dist_be <- dist_be[dist_be > -1]
  if (length(dist_be) > 0) {
    dist_be <- min(dist_be)
  } else {
    dist_be <- NA
  }
  min_dist_to_bound_end$ss[i] <- dist_be
}


for (j in 1:10000) {
  human_ss$random1_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
  dv <- c() # distance vector
  for (i in 1:nrow(human_ss)) {
    geneN <- human_ss[i,]$gene
    aa <- human_ss[i,]$aa
    r1 <- human_ss[i,]$random1_aa
    be <- dt[dt$gene == geneN]$boundary_end
    # r1
    r1_be <- r1 - be
    r1_be <- r1_be[r1_be > -1]
    if (length(r1_be) > 0) {
      r1_be <- min(r1_be)
    } else {
      r1_be <- NA
    }
    dv <- c(dv, r1_be)
  }
  min_dist_to_bound_end <- cbind(min_dist_to_bound_end, dv)
  print(j)
}

###################
# save(min_dist_to_bound_end, file = "../../DATA/domains/min_dist_to_bound_end.Rsave")
###################
### get 1000 columns averaged, calculate p-values
### PLOT
random1000 <- min_dist_to_bound_end[,2:1001]
variance_random <- apply(random1000, 1, function(x){var(x[complete.cases(x)])})
avg_random <- apply(random1000, 1, function(x){mean(x[complete.cases(x)])})

min_dist <- data.table(css = min_dist_to_bound_end$ss, control = avg_random)

data <- melt(min_dist)
#ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25) + xlim(0,300)
ggplot(data, aes(x=value, fill=variable)) + geom_density(alpha=0.25) + xlim(0,100)
ggplot(data, aes(x=value)) + geom_freqpoly(aes(colour = variable, linetype = variable)) + xlim(0,100) +
  scale_color_manual(values=c("#615388", "gray"), name = "position on mRNA", labels = c("stall site", "control")) +
  scale_linetype_manual(values = c("solid", "dotted")) + theme_classic()
#ggplot(data2, aes(x=value, fill=variable)) + geom_histogram(aes(color = variable), fill = "white", position = "dodge", binwidth = 1) + xlim(0,50)
# ggplot(data2, aes(x=value, fill=variable)) + stat_ecdf()

min_dist2 <- data.table(css = min_dist_to_bound_end$ss, random1 = min_dist_to_bound_end[,10], random2 = min_dist_to_bound_end[,100])
data2 <- melt(min_dist2)
ggplot(data2, aes(x=value)) + geom_freqpoly(aes(colour = variable, linetype = variable)) + xlim(0,100) +
  scale_color_manual(values=c("#615388", "black", "black"), name = "position on mRNA", labels = c("stall site", "random 1", "random2")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) + theme_classic()


############ subset?

# human_ss_subset <- human_ss[complete.cases(min_dist_to_bound_end),]
# random_pool_subset <- random_pool[names(random_pool) %in% human_ss_subset$h_tx]
# min_dist_to_bound_end <- min_dist_to_bound_end[complete.cases(min_dist_to_bound_end),]
#             
# for (j in 1:1000){
#   human_ss_subset$random1_aa <- sapply(as.character(human_ss_subset$h_tx), function(x){sample(random_pool_subset[[x]],1)})
#   distance_vector <- c()
#   for (i in 1:nrow(human_ss_subset)) {
#     geneN <- human_ss_subset[i,]$gene
#     aa <- human_ss_subset[i,]$aa
#     r1 <- human_ss_subset[i,]$random1_aa
#     be <- dt[dt$gene == geneN]$boundary_end
#     # r1
#     r1_be <- r1 - be
#     r1_be <- r1_be[r1_be > -1]
#     if (length(r1_be) > 0) {
#       r1_be <- min(r1_be)
#     } else {
#       r1_be <- NA
#     }
#     distance_vector <- c(distance_vector, r1_be)
#   }
#   min_dist_to_bound_end <- cbind(min_dist_to_bound_end, distance_vector)
#   print(j)
# }


# human_ss$random1_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# human_ss$random2_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# human_ss$random3_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# human_ss$random4_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# human_ss$random5_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# human_ss$random6_aa <- sapply(as.character(human_ss$h_tx), function(x){sample(random_pool[[x]],1)})
# 
# ###
# min_dist_to_bound_end <- data.table(gene = c(), dist = c(), r1 = c(), r2 = c(), r3 = c(), r4 = c(), r5 = c(), r6 = c())
# for (i in 1:nrow(human_ss)) {
#   geneN <- human_ss[i,]$gene
#   aa <- human_ss[i,]$aa
#   r1 <- human_ss[i,]$random1_aa
#   r2 <- human_ss[i,]$random2_aa
#   r3 <- human_ss[i,]$random3_aa
#   r4 <- human_ss[i,]$random4_aa
#   r5 <- human_ss[i,]$random5_aa
#   r6 <- human_ss[i,]$random6_aa
#   be <- dt[dt$gene == geneN]$boundary_end
#   # ss
#   dist_be <- aa - be
#   dist_be <- dist_be[dist_be > -1]
#   if (length(dist_be) > 0) {
#     dist_be <- min(dist_be)
#   } else {
#     dist_be <- NA
#   }
#   # r1
#   r1_be <- r1 - be
#   r1_be <- r1_be[r1_be > -1]
#   if (length(r1_be) > 0) {
#     r1_be <- min(r1_be)
#   } else {
#     r1_be <- NA
#   }
#   # r2
#   r2_be <- r2 - be
#   r2_be <- r2_be[r2_be > -1]
#   if (length(r2_be) > 0) {
#     r2_be <- min(r2_be)
#   } else {
#     r2_be <- NA
#   }
#   # r3
#   r3_be <- r3 - be
#   r3_be <- r3_be[r3_be > -1]
#   if (length(r3_be) > 0) {
#     r3_be <- min(r3_be)
#   } else {
#     r3_be <- NA
#   }
#   # r4
#   r4_be <- r4 - be
#   r4_be <- r4_be[r4_be > -1]
#   if (length(r4_be) > 0) {
#     r4_be <- min(r4_be)
#   } else {
#     r4_be <- NA
#   }
#   # r5
#   r5_be <- r5 - be
#   r5_be <- r5_be[r5_be > -1]
#   if (length(r5_be) > 0) {
#     r5_be <- min(r5_be)
#   } else {
#     r5_be <- NA
#   }
#   # r6
#   r6_be <- r6 - be
#   r6_be <- r6_be[r6_be > -1]
#   if (length(r6_be) > 0) {
#     r6_be <- min(r6_be)
#   } else {
#     r6_be <- NA
#   }
#   dist_be <- data.table(gene = geneN, dist = dist_be, r1 = r1_be, r2 = r2_be, r3 = r3_be, r4 = r4_be, r5 = r5_be, r6 = r6_be)
#   min_dist_to_bound_end <- rbind(min_dist_to_bound_end, dist_be)
#   print(i)
# }


data <- melt(min_dist_to_bound_end)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25) + xlim(0,300)

data2 <- data[data$variable == "dist" | data$variable == "r3"  | data$variable == "r1",]
ggplot(data2, aes(x=value, fill=variable)) + geom_density(alpha=0.25) + xlim(0,100)
#ggplot(data2, aes(x=value, fill=variable)) + geom_histogram(binwidth = 1, alpha=0.25) + xlim(0,100)
#ggplot(data2, aes(x=value, fill=variable)) + geom_area(aes(fill = variable), stat ="bin", alpha=0.6) + xlim(0,100)
ggplot(data2, aes(x=value)) + geom_freqpoly(aes(colour = variable, linetype = variable)) + xlim(0,100) +
  scale_color_manual(values=c("#615388", "gray", "gray"), name = "position on mRNA", labels = c("stall site", "random 1", "random2")) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted")) + theme_classic()
#ggplot(data2, aes(x=value, fill=variable)) + geom_histogram(aes(color = variable), fill = "white", position = "dodge", binwidth = 1) + xlim(0,50)
# ggplot(data2, aes(x=value, fill=variable)) + stat_ecdf()

data_dist <- data[data$variable == "dist",]
data_rand <- data[data$variable == "r3",]
ggplot(data_dist, aes(x=value, fill=variable)) + geom_histogram(binwidth = 1, alpha=1, fill = "#615388") +
  xlim(0,50) + ylim(0,35) + theme_bw()
ggplot(data_rand, aes(x=value, fill=variable)) + geom_histogram(binwidth = 1, alpha=1, fill = "gray") +
  xlim(0,50) + ylim(0,35) + theme_bw()

