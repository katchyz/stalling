### check GO_css and GO_control group:
# median ribo-seq count
# AA bias distribution
library(data.table)
library(ggplot2)
library(seqinr)

library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

########################################################################################################

load("../../DATA/stats/BIG_TABLE.Rsave")
setDT(BIG_TABLE)

# CSS
BIG_TABLE[css > 0]$ens_gene
# CTR: no peaks
BIG_TABLE[peaks == 0]$ens_gene
# CTR: non-conserved peaks
BIG_TABLE[peaks > 0 & css == 0]$ens_gene

# Ribo-seq
h1 <- get(load("../../DATA/peaks/human_h1.Rsave"))
h2 <- get(load("../../DATA/peaks/human_h2.Rsave"))
h3 <- get(load("../../DATA/peaks/human_h3.Rsave"))
h4 <- get(load("../../DATA/peaks/human_h4.Rsave"))

mean_ribo_h1 <- h1[, mean(riboseq), by = seqnames]
mean_ribo_h2 <- h2[, mean(riboseq), by = seqnames]
mean_ribo_h3 <- h3[, mean(riboseq), by = seqnames]
mean_ribo_h4 <- h4[, mean(riboseq), by = seqnames]

#mean_ribo_h1$fpkm <- (mean_ribo_h1$V1) * (10^9 / sum(h1$riboseq)) # not sure if sum is ok...
#mean_ribo_h2$fpkm <- (mean_ribo_h2$V1) * (10^9 / sum(h2$riboseq))
#mean_ribo_h3$fpkm <- (mean_ribo_h3$V1) * (10^9 / sum(h3$riboseq))
#mean_ribo_h4$fpkm <- (mean_ribo_h4$V1) * (10^9 / sum(h4$riboseq))

BIG_TABLE$ribo_h1 <- sapply(BIG_TABLE$tx, function(x){mean_ribo_h1[mean_ribo_h1$seqnames == x,]$V1})
BIG_TABLE$ribo_h2 <- sapply(BIG_TABLE$tx, function(x){mean_ribo_h2[mean_ribo_h2$seqnames == x,]$V1})
BIG_TABLE$ribo_h3 <- sapply(BIG_TABLE$tx, function(x){mean_ribo_h3[mean_ribo_h3$seqnames == x,]$V1})
BIG_TABLE$ribo_h4 <- sapply(BIG_TABLE$tx, function(x){mean_ribo_h4[mean_ribo_h4$seqnames == x,]$V1})

BIG_TABLE$ribo_h1 <- as.numeric(BIG_TABLE$ribo_h1)
BIG_TABLE$ribo_h2 <- as.numeric(BIG_TABLE$ribo_h2)
BIG_TABLE$ribo_h3 <- as.numeric(BIG_TABLE$ribo_h3)
BIG_TABLE$ribo_h4 <- as.numeric(BIG_TABLE$ribo_h4)

BIG_TABLE$set <- 'CSS'
BIG_TABLE[peaks == 0]$set <- 'CTR_no_peak'
BIG_TABLE[css == 0 & peaks > 0]$set <- 'CTR_not_conserved'

# plot
ggplot(BIG_TABLE, aes(x = ribo_h1, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))
# ggsave("mean_ribo_h1.png")
ggplot(BIG_TABLE, aes(x = ribo_h2, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))
# ggsave("mean_ribo_h2.png")
ggplot(BIG_TABLE, aes(x = ribo_h3, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))
# ggsave("mean_ribo_h3.png")
ggplot(BIG_TABLE, aes(x = ribo_h4, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))
# ggsave("mean_ribo_h4.png")

wilcox.test(BIG_TABLE[css > 0]$ribo_h1, BIG_TABLE[css == 0]$ribo_h1) # 1131 / 1412
wilcox.test(BIG_TABLE[css > 0]$ribo_h2, BIG_TABLE[css == 0]$ribo_h2) # 1514 / 1991
wilcox.test(BIG_TABLE[css > 0]$ribo_h3, BIG_TABLE[css == 0]$ribo_h3) # 1334 / 1650
wilcox.test(BIG_TABLE[css > 0]$ribo_h4, BIG_TABLE[css == 0]$ribo_h4) # 744 / 600

wilcox.test(BIG_TABLE[css > 0]$ribo_h1, BIG_TABLE[peaks == 0]$ribo_h1) # 1131 / 163
wilcox.test(BIG_TABLE[css > 0]$ribo_h2, BIG_TABLE[peaks == 0]$ribo_h2) # 1514 / 279
wilcox.test(BIG_TABLE[css > 0]$ribo_h3, BIG_TABLE[peaks == 0]$ribo_h3) # 1334 / 195
wilcox.test(BIG_TABLE[css > 0]$ribo_h4, BIG_TABLE[peaks == 0]$ribo_h4) # 744 / 67

wilcox.test(BIG_TABLE[css > 0]$ribo_h1, BIG_TABLE[peaks > 0 & css == 0]$ribo_h1) # 1131 / 1249
wilcox.test(BIG_TABLE[css > 0]$ribo_h2, BIG_TABLE[peaks > 0 & css == 0]$ribo_h2) # 1514 / 1712
wilcox.test(BIG_TABLE[css > 0]$ribo_h3, BIG_TABLE[peaks > 0 & css == 0]$ribo_h3) # 1334 / 1455
wilcox.test(BIG_TABLE[css > 0]$ribo_h4, BIG_TABLE[peaks > 0 & css == 0]$ribo_h4) # 744 / 533

l4_css <- c(BIG_TABLE[css > 0]$ribo_h1,
            BIG_TABLE[css > 0]$ribo_h2,
            BIG_TABLE[css > 0]$ribo_h3,
            BIG_TABLE[css > 0]$ribo_h4)

l4_ct1 <- c(BIG_TABLE[peaks == 0]$ribo_h1,
            BIG_TABLE[peaks == 0]$ribo_h2,
            BIG_TABLE[peaks == 0]$ribo_h3,
            BIG_TABLE[peaks == 0]$ribo_h4)

l4_ct2 <- c(BIG_TABLE[peaks > 0 & css == 0]$ribo_h1,
            BIG_TABLE[peaks > 0 & css == 0]$ribo_h2,
            BIG_TABLE[peaks > 0 & css == 0]$ribo_h3,
            BIG_TABLE[peaks > 0 & css == 0]$ribo_h4)

wilcox.test(l4_css, l4_ct1)
wilcox.test(l4_css, l4_ct2)
wilcox.test(l4_ct1, l4_ct2)

l4 <- data.table(ribo = c(l4_css[!is.na(l4_css)], l4_ct1[!is.na(l4_ct1)], l4_ct2[!is.na(l4_ct2)])) 
l4$set <- c(rep('CSS', length(l4_css[!is.na(l4_css)])),
            rep('control_no_peak', length(l4_ct1[!is.na(l4_ct1)])),
            rep('control_no_cons', length(l4_ct2[!is.na(l4_ct2)])))

#ggplot(l4, aes(x = ribo, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,10) + theme_classic() +
#  scale_fill_manual(values=c("#615388", "gray", "gray2"), name = "set",
#                    labels = c("CSS", "control: not conserved", "control: no peak"))
#ggsave("ribo_css_ct1_ct2.png")

# TE: Ribo-seq/RNA-seq
BIG_TABLE$rna_h2 <- sapply(BIG_TABLE$tx, function(x){cds_RNA_FPKM[names(cds_RNA_FPKM) == x]})

ggplot(BIG_TABLE, aes(x = rna_h2, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,100) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))

BIG_TABLE$ribo_fpkm_h2 <- BIG_TABLE$ribo_h2 * (10^9 / read_number_RIBO)

ggplot(BIG_TABLE, aes(x = ribo_rna_h2, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,50) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))
# ggsave("ribo_rna_h2.png")

wilcox.test(BIG_TABLE[css > 0]$ribo_fpkm_h2/BIG_TABLE[css > 0]$rna_h2,
            BIG_TABLE[peaks == 0]$ribo_fpkm_h2/BIG_TABLE[peaks == 0]$rna_h2)

wilcox.test(BIG_TABLE[css > 0]$ribo_fpkm_h2/BIG_TABLE[css > 0]$rna_h2,
            BIG_TABLE[peaks > 0 & css == 0]$ribo_fpkm_h2/BIG_TABLE[peaks > 0 & css == 0]$rna_h2)


### exclude lowly expressed from control, do GO again
wilcox.test(BIG_TABLE[css > 0 & !is.na(ribo_h2)]$ribo_h2, head(sort(BIG_TABLE[css == 0]$ribo_h2, decreasing = T), 1600))

bt <- BIG_TABLE[css == 0 & !is.na(ribo_h2)]
bt <- bt[order(-ribo_h2)][1:1600]

sim_ribo_MF <- enrichGO(gene     = bt$ens_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   universe      = BIG_TABLE$ens_gene,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

sim_ribo_BP <- enrichGO(gene     = bt$ens_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   universe      = BIG_TABLE$ens_gene,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(sim_ribo_MF, showCategory = 25, x = "Count")
dotplot(sim_ribo_BP, showCategory = 25, x = "Count")

BT_sim_ribo <- rbind(bt, BIG_TABLE[css > 0 & !is.na(ribo_h2)])
BT_sim_ribo$set2 <- c(rep("similar_ribo", nrow(bt)), rep("CSS", nrow(BIG_TABLE[css > 0 & !is.na(ribo_h2)])))
ggplot(BT_sim_ribo, aes(x = ribo_h2, fill = set2)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray"))
# ggsave("similar_ribo_h2.png")

### exclude low TE from control, do GO again
BIG_TABLE$ribo_rna_h2 <- BIG_TABLE$ribo_fpkm_h2 / BIG_TABLE$rna_h2
#wilcox.test(BIG_TABLE[css > 0]$ribo_rna_h2, head(sort(BIG_TABLE[css == 0]$ribo_rna_h2, decreasing = T), 1514))

#mean(BIG_TABLE[css > 0 & !is.na(ribo_rna_h2)]$ribo_rna_h2)
bte <- BIG_TABLE[css == 0 & !is.na(ribo_rna_h2)]
bte <- bte[order(ribo_rna_h2)][101:1550]
wilcox.test(BIG_TABLE[css > 0]$ribo_rna_h2, bte$ribo_rna_h2)

sim2_ribo_MF <- enrichGO(gene     = bte$ens_gene,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        universe      = BIG_TABLE$ens_gene,
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

sim2_ribo_BP <- enrichGO(gene     = bte$ens_gene,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        universe      = BIG_TABLE$ens_gene,
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

dotplot(sim2_ribo_MF, showCategory = 25, x = "Count")
dotplot(sim2_ribo_BP, showCategory = 25, x = "Count")

BT_sim_ribo_rna <- rbind(bte, BIG_TABLE[css > 0 & !is.na(ribo_rna_h2)])
BT_sim_ribo_rna$set2 <- c(rep("similar_TE", nrow(bte)), rep("CSS", nrow(BIG_TABLE[css > 0 & !is.na(ribo_rna_h2)])))
ggplot(BT_sim_ribo_rna, aes(x = ribo_rna_h2, fill = set2)) + geom_histogram(alpha=0.5) + xlim(0,20) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray"))
# ggsave("similar_TE_h2.png")



# AA distribution

fasta_css <- fasta[names(fasta) %in% BIG_TABLE[css > 0]$tx]
fasta_ct1 <- fasta[names(fasta) %in% BIG_TABLE[peaks == 0]$tx]
fasta_ct2 <- fasta[names(fasta) %in% BIG_TABLE[peaks > 0 & css == 0]$tx]

df <- data.frame(table(as.character(unlist(fasta_css))) / length(unlist(fasta_css)))
colnames(df) <- c('aa', 'freq_css')
df$freq_ct1 <- table(as.character(unlist(fasta_ct1))) / length(unlist(fasta_ct1))
df$freq_ct2 <- table(as.character(unlist(fasta_ct2))) / length(unlist(fasta_ct2))
freq_genome <- table(as.character(unlist(fasta))) / length(unlist(fasta))
df$freq_genome <- freq_genome[2:length(freq_genome)]
df <- df[-21,] #remove x
df <- df[-18,] #remove u
df$freq_css <- df$freq_css / sum(df$freq_css)
df$freq_ct1 <- df$freq_ct1 / sum(df$freq_ct1)
df$freq_ct2 <- df$freq_ct2 / sum(df$freq_ct2)
df$freq_gen <- df$freq_genome / sum(df$freq_genome)

chisq.test(df$freq_css*1000, p = df$freq_gen)
chisq.test(df$freq_ct1*1000, p = df$freq_gen)
chisq.test(df$freq_ct2*1000, p = df$freq_gen)

dt = data.table(aa = toupper(df$aa), freq = c(df$freq_css, df$freq_ct1, df$freq_ct2, df$freq_gen),
                type = c(rep("CSS", 20), rep("CTR_no_peak", 20), rep("CTR_not_conserved", 20), rep("genome", 20)))

ggplot(dt, aes(x = aa, y = freq, fill = type)) + geom_bar(stat = 'identity', position = 'dodge') + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray42", "black"))

# ggsave("AA_distribution.png", height = 6, width = 12)


