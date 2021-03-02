### alignment scores on CSSs vs non_conserved
library(data.table)
library(ggplot2)

load("../../DATA/stats/BIG_TABLE.Rsave")
setDT(BIG_TABLE)

###
css <- BIG_TABLE[css > 0]$gene_name # CSS # 1730
ctr_not_conserved <- BIG_TABLE[peaks > 0 & css == 0]$gene_name # CTR_not_conserved 2386

##################

cons_aln <- read.csv("../../DATA/GO/CONS_ALIGNMENTS.txt",
                     sep = '\t', header = F)
colnames(cons_aln) <- c("filename", "org1", "org2", "leven")
setDT(cons_aln)
cons_aln$gene <- sapply(as.character(cons_aln$filename), function(x){unlist(strsplit(x, split = "_"))[1]})

sum(cons_aln$gene %in% css) # 5363
sum(cons_aln$gene %in% ctr_not_conserved) # 4418

# get human only
cons_aln$org1_org <- sapply(as.character(cons_aln$org1), function(x){unlist(strsplit(x, split = "_"))[1]})
cons_aln$org2_org <- sapply(as.character(cons_aln$org2), function(x){unlist(strsplit(x, split = "_"))[1]})

cons_aln_h <- cons_aln[org1_org == "human" | org2_org == "human"]
cons_aln_hm <- cons_aln[(org1_org == "mouse" & org2_org == "human") | (org1_org == "human" & org2_org == "mouse")]
cons_aln_hm <- cons_aln_hm[, .SD[1], gene]

# get max
#css_max <- cons_aln_hm[gene %in% css]
#ctr_max <- cons_aln_hm[gene %in% ctr_not_conserved]
#
#dt_css <- css_max[, max(leven), by = gene]
#colnames(dt_css) <- c("gene", "leven")
#dt_css$set <- "CSS"
#
#dt_ctr <- ctr_max[, max(leven), by = gene]
#colnames(dt_ctr) <- c("gene", "leven")
#dt_ctr$set <- "CTR_not_conserved"

dt_css <- data.table(gene = cons_aln_h[gene %in% css]$gene, leven = cons_aln_h[gene %in% css]$leven, set = 'CSS')
dt_ctr <- data.table(gene = cons_aln_h[gene %in% ctr_not_conserved]$gene,
                     leven = cons_aln_h[gene %in% ctr_not_conserved]$leven, set = 'CTR_not_conserved')
dt <- rbind(dt_css, dt_ctr)

dt2 <- dt[!(set == "CTR_not_conserved" & leven < 64.2)]

ggplot(dt2, aes(x = leven, fill = set)) + geom_histogram(alpha=0.5) + xlim(50,100) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray"), labels=c("CSS", "similar_aln"))
ggsave("../../DATA/GO/similar_aln_score.png")

wilcox.test(dt2[set == "CSS"]$leven, dt2[set == "CTR_not_conserved"]$leven)

median(dt[set == "CSS"]$leven) # 82.8%
median(dt[set == "CTR_not_conserved"]$leven) # 81.6%

# sample so that distributions are the same, do GO again
wilcox.test(dt[set == "CSS"]$leven, dt[set == "CTR_not_conserved" & leven > 64.2]$leven) # 0.99

ctr_aln <- dt[set == "CTR_not_conserved" & leven > 64.2]$gene
ctr_aln <- BIG_TABLE[gene_name %in% ctr_aln]$ens_gene

library(clusterProfiler)
library(org.Hs.eg.db)
go_path <- "../../DATA/GO"
load(file = c(file.path(go_path, "human_well_expressed_homologs.Rsave"))) # wet_homo

aln_MF <- enrichGO(gene          = ctr_aln,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wet_homo,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

aln_BP <- enrichGO(gene          = ctr_aln,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wet_homo,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(aln_MF, showCategory = 35, x = "Count")
dotplot(aln_BP, showCategory = 35, x = "Count")


### check for number of peaks
BIG_TABLE$set <- "CTR_no_peak"
BIG_TABLE[css > 0]$set <- "CSS"
BIG_TABLE[peaks > 0 & css == 0]$set <- "CTR_not_conserved"

ggplot(BIG_TABLE, aes(x = peaks, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,50) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "gray2"))

wilcox.test(BIG_TABLE[set == "CSS"]$peaks, BIG_TABLE[set == "CTR_not_conserved"]$peaks)
mean(BIG_TABLE[set == "CSS"]$peaks) # 13.9
mean(BIG_TABLE[set == "CTR_not_conserved"]$peaks) # 6.7

median(BIG_TABLE[set == "CTR_not_conserved" & peaks > 4]$peaks)
to_excl <- c(sample(BIG_TABLE[set == "CTR_not_conserved" & peaks == 4]$ens_gene, 70),
             sample(BIG_TABLE[set == "CTR_not_conserved" & peaks == 5]$ens_gene, 70))
wilcox.test(BIG_TABLE[set == "CSS"]$peaks,
            BIG_TABLE[set == "CTR_not_conserved" & peaks > 3 & !(ens_gene %in% to_excl)]$peaks) # 0.92


dt3 <- BIG_TABLE[set == "CSS"]
dt3 <- rbind(dt3, BIG_TABLE[set == "CTR_not_conserved" & peaks > 3 & !(ens_gene %in% to_excl)])

ggplot(dt3, aes(x = peaks, fill = set)) + geom_histogram(alpha=0.5) + xlim(0,50) + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray"), labels=c("CSS", "similar_peak_no"))
ggsave("../../DATA/GO/similar_peak_number.png")

ctr_peak <- BIG_TABLE[set == "CTR_not_conserved" & peaks > 3 & !(ens_gene %in% to_excl)]$ens_gene

peak_MF <- enrichGO(gene          = ctr_peak,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   universe      = wet_homo,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

peak_BP <- enrichGO(gene          = ctr_peak,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   universe      = wet_homo,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(peak_MF, showCategory = 35, x = "Count")
dotplot(peak_BP, showCategory = 35, x = "Count")


