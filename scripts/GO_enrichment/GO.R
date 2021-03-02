## GO analysis

# OrgDb
library(clusterProfiler)
library(org.Sc.sgd.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

library(ggplot2)
theme_set(theme_bw())
library(reshape2)
library(RColorBrewer)


#############################
#############################
## all homologs
## homologs with stall sites
## homologs without stall sites
library(GenomicFeatures)
library(plyr)
library(data.table)

gtf_path <- "../../DATA/genomes/GTF"

#### homologs with stall sites
css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")

txdb_y <- makeTxDbFromGFF(c(file.path(gtf_path, "Saccharomyces_cerevisiae.R64-1-1.79.gtf")), format="gtf")
txLengths_y <- transcriptLengths(txdb_y, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_y <- txLengths_y[order(txLengths_y$tx_name),]

txdb_f <- makeTxDbFromGFF(c(file.path(gtf_path, "Drosophila_melanogaster.BDGP6.79.gtf")), format="gtf")
txLengths_f <- transcriptLengths(txdb_f, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_f <- txLengths_f[order(txLengths_f$tx_name),]

txdb_z <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.81_chr.gtf")), format="gtf")
txLengths_z <- transcriptLengths(txdb_z, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_z <- txLengths_z[order(txLengths_z$tx_name),]

txdb_m <- makeTxDbFromGFF(c(file.path(gtf_path, "Mus_musculus.GRCm38.79.chr.gtf")), format="gtf")
txLengths_m <- transcriptLengths(txdb_m, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_m <- txLengths_m[order(txLengths_m$tx_name),]

txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]


ss_y <- as.character(css[!is.na(css$y_ss),]$y_tx) # stall sites
ss_y <- unique(ss_y) # genes
ss_y <- as.character(sapply(ss_y, function(x){txLengths_y[txLengths_y$tx_name == x,]$gene_id}))
  
ss_f <- as.character(css[!is.na(css$f_ss),]$f_tx) # stall sites
ss_f <- unique(ss_f) # genes
ss_f <- as.character(sapply(ss_f, function(x){txLengths_f[txLengths_f$tx_name == x,]$gene_id}))

ss_z <- as.character(css[!is.na(css$z_ss),]$z_tx) # stall sites
ss_z <- unique(ss_z) # genes
ss_z <- as.character(sapply(ss_z, function(x){txLengths_z[txLengths_z$tx_name == x,]$gene_id}))

ss_m <- as.character(css[!is.na(css$m_ss),]$m_tx) # stall sites
ss_m <- unique(ss_m) # genes
ss_m <- as.character(sapply(ss_m, function(x){txLengths_m[txLengths_m$tx_name == x,]$gene_id}))

ss_h <- as.character(css[!is.na(css$h_ss),]$h_tx) # stall sites
ss_h <- unique(ss_h) # genes
ss_h <- as.character(sapply(ss_h, function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id}))



# yeast 
ss_Y_MF <- enrichGO(gene            = ss_y,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = 'ENSEMBL', # 'ORF'
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_Y_BP <- enrichGO(gene          = ss_y,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = 'ENSEMBL', # 'ORF'
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# fruit fly
ss_F_MF <- enrichGO(gene          = ss_f,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_F_BP <- enrichGO(gene          = ss_f,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(ss_F_MF, showCategory = 15)
dotplot(ss_F_BP, showCategory = 15)

# zebrafish
ss_Z_MF <- enrichGO(gene          = ss_z,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_Z_BP <- enrichGO(gene          = ss_z,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(ss_Z_MF, showCategory = 15)
dotplot(ss_Z_BP, showCategory = 15)

# mouse
ss_M_MF <- enrichGO(gene          = ss_m,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_M_BP <- enrichGO(gene          = ss_m,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(ss_M_MF, showCategory = 15)
dotplot(ss_M_BP, showCategory = 15)

# human
ss_H_MF <- enrichGO(gene          = ss_h,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_H_BP <- enrichGO(gene          = ss_h,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ss_H_CC <- enrichGO(gene          = ss_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(ss_H_MF, showCategory = 15)
dotplot(ss_H_BP, showCategory = 15)
dotplot(ss_H_CC, showCategory = 15)

#### homologs without stall sites
no_y <- homo_y[!homo_y$ens_gene %in% ss_y,]$ens_gene
no_f <- homo_f[!homo_f$ens_gene %in% ss_f,]$ens_gene
no_z <- homo_z[!homo_z$ens_gene %in% ss_z,]$ens_gene
no_m <- homo_m[!homo_m$ens_gene %in% ss_m,]$ens_gene
no_h <- homo_h[!homo_h$ens_gene %in% ss_h,]$ens_gene


# yeast 
no_Y_MF <- enrichGO(gene            = no_y,
                    OrgDb         = org.Sc.sgd.db,
                    keyType       = 'ENSEMBL', # 'ORF'
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_Y_BP <- enrichGO(gene          = no_y,
                    OrgDb         = org.Sc.sgd.db,
                    keyType       = 'ENSEMBL', # 'ORF'
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

# fruit fly
no_F_MF <- enrichGO(gene          = no_f,
                    OrgDb         = org.Dm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_F_BP <- enrichGO(gene          = no_f,
                    OrgDb         = org.Dm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(no_F_MF, showCategory = 15)
dotplot(no_F_BP, showCategory = 15)

# zebrafish
no_Z_MF <- enrichGO(gene          = no_z,
                    OrgDb         = org.Dr.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_Z_BP <- enrichGO(gene          = no_z,
                    OrgDb         = org.Dr.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(no_Z_MF, showCategory = 15)
dotplot(no_Z_BP, showCategory = 15)

# mouse
no_M_MF <- enrichGO(gene          = no_m,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_M_BP <- enrichGO(gene          = no_m,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(no_M_MF, showCategory = 15)
dotplot(no_M_BP, showCategory = 15)

# human
no_H_MF <- enrichGO(gene          = no_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_H_BP <- enrichGO(gene          = no_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(ss_H_MF, showCategory = 15)
dotplot(ss_H_BP, showCategory = 15)


####
ss_H_MF$Description
no_H_MF$Description
homo_H_MF$Description

ss_only_H_MF <- ss_H_MF$Description[!ss_H_MF$Description %in% homo_H_MF$Description]
ss_only_H_MF <- ss_only_H_MF[!ss_only_H_MF %in% no_H_MF$Description]

ss_H_MF[ss_H_MF$Description %in% ss_only_H_MF]$GeneRatio





# background
# yeast
bckg_Y_MF <- enrichGO(gene          = ss_y,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = homo_y$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_Y_BP <- enrichGO(gene          = ss_y,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = homo_y$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)


# fruitfly
bckg_F_MF <- enrichGO(gene          = ss_f,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = homo_f$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_F_BP <- enrichGO(gene          = ss_f,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = homo_f$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)


# zebrafish
bckg_Z_MF <- enrichGO(gene          = ss_z,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = homo_z$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_Z_BP <- enrichGO(gene          = ss_z,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = homo_z$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)


# mouse
bckg_M_MF <- enrichGO(gene          = ss_m,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = homo_m$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_M_BP <- enrichGO(gene          = ss_m,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = homo_m$ens_gene,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# human
bckg_H_MF <- enrichGO(gene          = ss_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = unique(homo_h$ens_gene),
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_H_BP <- enrichGO(gene          = ss_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = unique(homo_h$ens_gene),
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)


dotplot(bckg_F_MF, showCategory = 15)
dotplot(bckg_F_BP, showCategory = 15)
dotplot(bckg_Z_MF, showCategory = 15)
dotplot(bckg_Z_BP, showCategory = 15)
dotplot(bckg_M_MF, showCategory = 15)
dotplot(bckg_M_BP, showCategory = 15)
dotplot(bckg_H_MF, showCategory = 15)
dotplot(bckg_H_BP, showCategory = 15)


# background - no stall sites
# yeast 
bckg_no_Y_MF <- enrichGO(gene            = no_y,
                    OrgDb         = org.Sc.sgd.db,
                    keyType       = 'ENSEMBL', # 'ORF'
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = homo_y$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_no_Y_BP <- enrichGO(gene          = no_y,
                    OrgDb         = org.Sc.sgd.db,
                    keyType       = 'ENSEMBL', # 'ORF'
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = homo_y$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

# fruit fly
bckg_no_F_MF <- enrichGO(gene          = no_f,
                    OrgDb         = org.Dm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = homo_f$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_no_F_BP <- enrichGO(gene          = no_f,
                    OrgDb         = org.Dm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = homo_f$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(bckg_no_F_MF, showCategory = 15)
dotplot(bckg_no_F_BP, showCategory = 15)

# zebrafish
bckg_no_Z_MF <- enrichGO(gene          = no_z,
                    OrgDb         = org.Dr.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = homo_z$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_no_Z_BP <- enrichGO(gene          = no_z,
                    OrgDb         = org.Dr.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = homo_z$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(bckg_no_Z_MF, showCategory = 15)
dotplot(bckg_no_Z_BP, showCategory = 15)

# mouse
bckg_no_M_MF <- enrichGO(gene          = no_m,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = homo_m$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_no_M_BP <- enrichGO(gene          = no_m,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = homo_m$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(bckg_no_M_MF, showCategory = 15)
dotplot(bckg_no_M_BP, showCategory = 15)

# human
bckg_no_H_MF <- enrichGO(gene          = no_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    universe      = homo_h$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

bckg_no_H_BP <- enrichGO(gene          = no_h,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = homo_h$ens_gene,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(bckg_ss_H_MF, showCategory = 15)
dotplot(bckg_ss_H_BP, showCategory = 15)

############################################

go_mf <- c("ss_Y_MF", "ss_F_MF", "ss_Z_MF", "ss_M_MF", "ss_H_MF", # stall site against genome
           "no_Y_MF", "no_F_MF", "no_Z_MF", "no_M_MF", "no_H_MF", # no stall site against genome
           "homo_Y_MF", "homo_F_MF", "homo_Z_MF", "homo_M_MF", "homo_H_MF", # homologs against genome
           "bckg_Y_MF", "bckg_F_MF", "bckg_Z_MF", "bckg_M_MF", "bckg_H_MF") # stall sites against homologs

go_bp <- c("ss_Y_BP", "ss_F_BP", "ss_Z_BP", "ss_M_BP", "ss_H_BP", # stall site against genome
           "no_Y_BP", "no_F_BP", "no_Z_BP", "no_M_BP", "no_H_BP", # no stall site against genome
           "homo_Y_BP", "homo_F_BP", "homo_Z_BP", "homo_M_BP", "homo_H_BP", # homologs against genome
           "bckg_Y_BP", "bckg_F_BP", "bckg_Z_BP", "bckg_M_BP", "bckg_H_BP") # stall sites against homologs


### just human
dotplot(ss_H_MF, showCategory = 30) # genes: 1674, bckg: 19626
dotplot(no_H_MF, showCategory = 30) # genes: 14367, bckg: 19626
dotplot(homo_H_MF, showCategory = 30) # genes: 16040, bckg: 19626
dotplot(bckg_H_MF, showCategory = 30) # genes: 1674, bckg: 16040

dotplot(ss_H_BP, showCategory = 30)
dotplot(no_H_BP, showCategory = 30)
dotplot(homo_H_BP, showCategory = 30)
dotplot(bckg_H_BP, showCategory = 30)

# check well-expressed genes and expression of other homologs

wet <- get(load(file = "../../DATA/peaks/human_h1.Rsave"))
wet2 <- get(load(file = "../../DATA/peaks/human_h2.Rsave"))
wet3 <- get(load(file = "../../DATA/peaks/human_h3.Rsave")) # has "ENST00000368885"
wet4 <- get(load(file = "../../DATA/peaks/human_h4.Rsave")) # has "ENST00000368885"

wet_tx <- unique(as.character(wet$seqnames))
wet_gene <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx,]$gene_id)
wet_tx2 <- unique(as.character(wet2$seqnames))
wet_gene2 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx2,]$gene_id)
wet_tx3 <- unique(as.character(wet3$seqnames))
wet_gene3 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx3,]$gene_id)
wet_tx4 <- unique(as.character(wet4$seqnames))
wet_gene4 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx4,]$gene_id)

wet_t <- rbind(data.table(seqnames = wet$seqnames, zscore = wet$zscore),
               data.table(seqnames = wet2$seqnames, zscore = wet2$zscore))
wet_t <- rbind(wet_t,
               data.table(seqnames = wet3$seqnames, zscore = wet3$zscore))
wet_t <- rbind(wet_t,
               data.table(seqnames = wet4$seqnames, zscore = wet4$zscore))

wet_g <- c(wet_gene, wet_gene2, wet_gene3, wet_gene4)
wet_g <- unique(wet_g)

ss_h_unique <- unique(ss_h)
no_h_unique <- unique(no_h)

ss_h_unique <- ss_h_unique[ss_h_unique %in% wet_g]
no_h_unique <- no_h_unique[no_h_unique %in% wet_g]

# human
ss_wet_H_MF <- enrichGO(gene          = ss_h_unique,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

ss_wet_H_BP <- enrichGO(gene          = ss_h_unique,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)


no_wet_H_MF <- enrichGO(gene          = no_h_unique,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_wet_H_BP <- enrichGO(gene          = no_h_unique,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

dotplot(ss_wet_H_MF, showCategory = 30)
dotplot(ss_wet_H_BP, showCategory = 30)
dotplot(no_wet_H_MF, showCategory = 30)
dotplot(no_wet_H_BP, showCategory = 30)

wet_H_MF <- enrichGO(gene          = wet_g,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

wet_H_BP <- enrichGO(gene          = wet_g,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

dotplot(wet_H_MF, showCategory = 30)
dotplot(wet_H_BP, showCategory = 30)


##########################################################
### export list ss_H_MF and SS_H_BP
### compare to no_wet_H_MF and no_wet_H_BP

# ss_H_MF$Description[!ss_H_MF$Description %in% lowz_H_MF$Description]
# ss_H_BP$Description[!ss_H_BP$Description %in% lowz_H_BP$Description]
# 
# dt_MF_human <- data.table(GO = ss_H_MF$Description, genes = ss_H_MF$geneID)
# dt_BP_human <- data.table(GO = ss_H_BP$Description, genes = ss_H_BP$geneID)
# 
# write.table(dt_MF_human, "../../DATA/GO/GO_peaks_human_MF.txt", row.names = F)
# write.table(dt_BP_human, "../../DATA/GO/GO_peaks_human_BP.txt", row.names = F)

### export BED file with human conserved stall sites
# library(GenomicRanges)
# library(GenomicFeatures)
# library(rtracklayer)
# 
# human_ss <- data.table(tx = css[!is.na(css$h_ss),]$h_tx, ss = css[!is.na(css$h_ss),]$h_ss)
# human_ss$end <- human_ss$ss + 2
# 
# gr <- makeGRangesFromDataFrame(human_ss,
#                          keep.extra.columns=FALSE,
#                          ignore.strand=FALSE,
#                          seqinfo=NULL,
#                          seqnames.field="tx",
#                          start.field="ss",
#                          end.field="end",
#                          starts.in.df.are.0based=FALSE)
# 
# txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
# exons_h <- exonsBy(txdb_h, by="tx", use.names=TRUE)
# 
# ss_genomic <- mapFromTranscripts(gr, exons_h)

# export.bed(ss_genomic, "../../DATA/GO/GO_peaks_human.bed") ### wrong coord.! nt in cds


