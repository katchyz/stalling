## GO analysis

# OrgDb
library(clusterProfiler)
library(org.Sc.sgd.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# genes with stall sites
go_path <- "/Volumes/USELESS/STALLING/GO"
yeast_ss <- as.character(read.table(file = c(file.path(go_path, "yeast_ss.txt")), header = F)$V1)
fruitfly_ss <- as.character(read.table(file = c(file.path(go_path, "fruitfly_ss.txt")), header = F)$V1)
zebrafish_ss <- as.character(read.table(file = c(file.path(go_path, "zebrafish_ss.txt")), header = F)$V1)
mouse_ss <- as.character(read.table(file = c(file.path(go_path, "mouse_ss.txt")), header = F)$V1)
human_ss <- as.character(read.table(file = c(file.path(go_path, "human_ss.txt")), header = F)$V1)


# over-representation analysis

Y_MF <- enrichGO(gene          = yeast_ss,
                 OrgDb         = org.Sc.sgd.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

Y_BP <- enrichGO(gene          = yeast_ss,
                 OrgDb         = org.Sc.sgd.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

F_MF <- enrichGO(gene          = fruitfly_ss,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

F_BP <- enrichGO(gene          = fruitfly_ss,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

Z_MF <- enrichGO(gene          = zebrafish_ss,
                 OrgDb         = org.Dr.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

Z_BP <- enrichGO(gene          = zebrafish_ss,
                 OrgDb         = org.Dr.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

M_MF <- enrichGO(gene          = mouse_ss,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

M_BP <- enrichGO(gene          = mouse_ss,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

H_MF <- enrichGO(gene          = human_ss,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

H_BP <- enrichGO(gene          = human_ss,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

dotplot(Y_MF, showCategory = 15) # 17
dotplot(F_MF, showCategory = 15) # 15
dotplot(Z_MF, showCategory = 15) # 0
dotplot(M_MF, showCategory = 15) # 41
dotplot(H_MF, showCategory = 15) # 31

dotplot(Y_BP, showCategory = 15) # 102
dotplot(F_BP, showCategory = 15) # 374
dotplot(Z_BP, showCategory = 15) # 19
dotplot(M_BP, showCategory = 15) # 207
dotplot(H_BP, showCategory = 15) # 223

# y.df <- bitr(yeast_ss, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Sc.sgd.db)
# f.df <- bitr(fruitfly_ss, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
# z.df <- bitr(zebrafish_ss, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
# m.df <- bitr(mouse_ss, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# h.df <- bitr(human_ss, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 
# cluster <- list(yeast = y.df$ENTREZID, fruitfly = f.df$ENTREZID, zebrafish = z.df$ENTREZID,
#                 mouse = m.df$ENTREZID, human = h.df$ENTREZID)
# 
# ck <- compareCluster(geneCluster = cluster, fun = "enrichGO")

# molecular function
mf <- c(Y_MF$ID, F_MF$ID, Z_MF$ID, M_MF$ID, H_MF$ID)
sort(table(mf), decreasing = T)

mf_desc <- c(Y_MF$Description, F_MF$Description, Z_MF$Description, M_MF$Description, H_MF$Description)
GOmf <- sort(table(mf_desc), decreasing = T)
GOmf[GOmf > 2]

# biological process
bp <- c(Y_BP$ID, F_BP$ID, Z_BP$ID, M_BP$ID, H_BP$ID)
sort(table(bp), decreasing = T)

bp_desc <- c(Y_BP$Description, F_BP$Description, Z_BP$Description, M_BP$Description, H_BP$Description)
GObp <- sort(table(bp_desc), decreasing = T)
GObp[GObp > 3]

# make bar/dot plot of common categories
# (with the proportion of genes in total of 5 organisms)

library(ggplot2)
theme_set(theme_bw())
library(reshape2)
library(RColorBrewer)

##### mf
mf <- data.frame(go = names(GOmf[GOmf > 2]))

get_GeneRatio <- function(go_term, enrichResult) {
  if (go_term %in% enrichResult$Description) {
    idx <- which(enrichResult$Description %in% go_term)
    gr <- enrichResult$GeneRatio[idx]
    gr <- as.numeric(strsplit(gr, split = "/")[[1]][1])
    return(gr)
  } else {
    return(0)
  }
}

ymf <- c()
for (go_term in mf$go) {
  gr <- get_GeneRatio(go_term, Y_MF)
  ymf <- c(ymf, gr)
}
mf$yeast <- ymf

fmf <- c()
for (go_term in mf$go) {
  gr <- get_GeneRatio(go_term, F_MF)
  fmf <- c(fmf, gr)
}
mf$fruitfly <- fmf

zmf <- c()
for (go_term in mf$go) {
  gr <- get_GeneRatio(go_term, Z_MF)
  zmf <- c(zmf, gr)
}
mf$zebrafish <- zmf

mmf <- c()
for (go_term in mf$go) {
  gr <- get_GeneRatio(go_term, M_MF)
  mmf <- c(mmf, gr)
}
mf$mouse <- mmf

hmf <- c()
for (go_term in mf$go) {
  gr <- get_GeneRatio(go_term, H_MF)
  hmf <- c(hmf, gr)
}
mf$human <- hmf

###### bp
bp <- data.frame(go = names(GObp[GObp > 3]))

ybp <- c()
for (go_term in bp$go) {
  gr <- get_GeneRatio(go_term, Y_BP)
  ybp <- c(ybp, gr)
}
bp$yeast <- ybp

fbp <- c()
for (go_term in bp$go) {
  gr <- get_GeneRatio(go_term, F_BP)
  fbp <- c(fbp, gr)
}
bp$fruitfly <- fbp

zbp <- c()
for (go_term in bp$go) {
  gr <- get_GeneRatio(go_term, Z_BP)
  zbp <- c(zbp, gr)
}
bp$zebrafish <- zbp

mbp <- c()
for (go_term in bp$go) {
  gr <- get_GeneRatio(go_term, M_BP)
  mbp <- c(mbp, gr)
}
bp$mouse <- mbp

hbp <- c()
for (go_term in bp$go) {
  gr <- get_GeneRatio(go_term, H_BP)
  hbp <- c(hbp, gr)
}
bp$human <- hbp

####### plot mf and bp as bar plots (diff org diff color in a bar)
# mf: yeast - 646, fruitfly - 690, zebrafish - 147, mouse - 479, human - 509
# total mf: 2471
# bp: yeast - 646, fruitfly - 711, zebrafish - 108, mouse - 486, human - 513
# toatl bp: 2464

mf_melt <- melt(mf)
ggplot(mf_melt, aes(x = go, y = value, fill = variable)) + geom_bar(stat = "identity") + coord_flip() +
  scale_fill_brewer(palette="Accent") + labs(fill = "organism") + theme(axis.title.y = element_blank()) +
  ylab("number of genes") + theme(axis.text.y = element_text(face = "bold", size = 15))

bp_temp <- bp
bp_temp$sum <- bp$yeast + bp$fruitfly + bp$zebrafish + bp$mouse + bp$human
bp_temp <- bp_temp[-5,]
bp_temp <- bp_temp[bp_temp$sum > 99,]
bp_temp$sum <- NULL
bp_melt <- melt(bp_temp)
ggplot(bp_melt, aes(x = go, y = value, fill = variable)) + geom_bar(stat = "identity") + coord_flip() +
  scale_fill_brewer(palette="Accent") + labs(fill = "organism") + theme(axis.title.y = element_blank()) +
  ylab("number of genes") + theme(axis.text.y = element_text(face = "bold", size = 12))


#############################
#############################
## all homologs
## homologs with stall sites
## homologs without stall sites
library(GenomicFeatures)
library(plyr)
library(data.table)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"

homologs_path <- "/Volumes/USELESS/STALLING/conservation/biomart_export"
#yeast
homo <- as.character(read.table(file = c(file.path(homologs_path, "yeast.txt")))$V1)
homo_tx <- as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][1]}))
homo_gn <- toupper(as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][2]})))
yeast <- data.table(tx_yeast = homo_tx, gene = homo_gn)
yeast <- yeast[!(yeast$gene == "N/A")]
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Saccharomyces_cerevisiae.R64-1-1.79.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
yeast <- yeast[order(yeast$tx_yeast),]
yeast <- yeast[yeast$tx_yeast %in% txLengths$tx_name]
txLengths <- txLengths[txLengths$tx_name %in% yeast$tx_yeast,]
yeast$ens_gene <- txLengths$gene_id
colnames(yeast) <- c("tx_yeast", "gene", "ens_gene_yeast")

#fruitfly
homo <- as.character(read.table(file = c(file.path(homologs_path, "fruitfly.txt")))$V1)
homo_tx <- as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][1]}))
homo_gn <- toupper(as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][2]})))
fruitfly <- data.table(tx_fruitfly = homo_tx, gene = homo_gn)
fruitfly <- fruitfly[!(fruitfly$gene == "N/A")]
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Drosophila_melanogaster.BDGP6.79.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
fruitfly <- fruitfly[order(fruitfly$tx_fruitfly),]
fruitfly <- fruitfly[fruitfly$tx_fruitfly %in% txLengths$tx_name]
txLengths <- txLengths[txLengths$tx_name %in% fruitfly$tx_fruitfly,]
fruitfly$ens_gene <- txLengths$gene_id
fruitfly$tx_len <- txLengths$tx_len
fruitfly <- arrange(fruitfly, ens_gene, desc(tx_len))
fruitfly <- fruitfly[!duplicated(fruitfly$ens_gene),]
colnames(fruitfly) <- c("tx_fruitfly", "gene", "ens_gene_fruitfly", "tx_len_f")

#zebrafish
homo <- as.character(read.table(file = c(file.path(homologs_path, "zebrafish.txt")))$V1)
homo_tx <- as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][1]}))
homo_gn <- toupper(as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][2]})))
zebrafish <- data.table(tx_zebrafish = homo_tx, gene = homo_gn)
zebrafish <- zebrafish[!(zebrafish$gene == "N/A")]
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.81_chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
zebrafish <- zebrafish[order(zebrafish$tx_zebrafish),]
zebrafish <- zebrafish[zebrafish$tx_zebrafish %in% txLengths$tx_name]
txLengths <- txLengths[txLengths$tx_name %in% zebrafish$tx_zebrafish,]
zebrafish$ens_gene <- txLengths$gene_id
zebrafish$tx_len <- txLengths$tx_len
zebrafish <- arrange(zebrafish, ens_gene, desc(tx_len))
zebrafish <- zebrafish[!duplicated(zebrafish$ens_gene),]
colnames(zebrafish) <- c("tx_zebrafish", "gene", "ens_gene_zebrafish", "tx_len_z")

#mouse
homo <- as.character(read.table(file = c(file.path(homologs_path, "mouse.txt")))$V1)
homo_tx <- as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][1]}))
homo_gn <- toupper(as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][2]})))
mouse <- data.table(tx_mouse = homo_tx, gene = homo_gn)
mouse <- mouse[!(mouse$gene == "N/A")]
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Mus_musculus.GRCm38.79.chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
mouse <- mouse[order(mouse$tx_mouse),]
mouse <- mouse[mouse$tx_mouse %in% txLengths$tx_name]
txLengths <- txLengths[txLengths$tx_name %in% mouse$tx_mouse,]
mouse$ens_gene <- txLengths$gene_id
mouse$tx_len <- txLengths$tx_len
mouse <- arrange(mouse, ens_gene, desc(tx_len))
mouse <- mouse[!duplicated(mouse$ens_gene),]
colnames(mouse) <- c("tx_mouse", "gene", "ens_gene_mouse", "tx_len_m")

#human
homo <- as.character(read.table(file = c(file.path(homologs_path, "human.txt")))$V1)
homo_tx <- as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][1]}))
homo_gn <- toupper(as.character(sapply(homo, function(x){strsplit(x, split = ",")[[1]][2]})))
human <- data.table(tx_human = homo_tx, gene = homo_gn)
human <- human[!(human$gene == "N/A")]
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
human <- human[order(human$tx_human),]
human <- human[human$tx_human %in% txLengths$tx_name]
txLengths <- txLengths[txLengths$tx_name %in% human$tx_human,]
human$ens_gene <- txLengths$gene_id
human$tx_len <- txLengths$tx_len
human <- arrange(human, ens_gene, desc(tx_len))
human <- human[!duplicated(human$ens_gene),]
colnames(human) <- c("tx_human", "gene", "ens_gene_human", "tx_len_h")

##

dt_all <- merge(yeast[complete.cases(yeast)], fruitfly[complete.cases(fruitfly)], by="gene", all=TRUE)
dt_all <- merge(dt_all, zebrafish[complete.cases(zebrafish)], by="gene", all=TRUE)
dt_all <- merge(dt_all, mouse[complete.cases(mouse)], by="gene", all=TRUE)
dt_all <- merge(dt_all, human[complete.cases(human)], by="gene", all=TRUE)

homo_y <- dt_all[(!is.na(dt_all$ens_gene_yeast)) & ( (!is.na(dt_all$ens_gene_fruitfly)) |
                                                       (!is.na(dt_all$ens_gene_zebrafish)) |
                                                          (!is.na(dt_all$ens_gene_mouse)) |
                                                             (!is.na(dt_all$ens_gene_human)) )]
homo_y <- data.table(gene = homo_y$gene, ens_gene = homo_y$ens_gene_yeast)


homo_f <- dt_all[(!is.na(dt_all$ens_gene_fruitfly)) & ( (!is.na(dt_all$ens_gene_yeast)) |
                                                       (!is.na(dt_all$ens_gene_zebrafish)) |
                                                       (!is.na(dt_all$ens_gene_mouse)) |
                                                       (!is.na(dt_all$ens_gene_human)) )]
homo_f <- data.table(gene = homo_f$gene, ens_gene = homo_f$ens_gene_fruitfly)

homo_z <- dt_all[(!is.na(dt_all$ens_gene_zebrafish)) & ( (!is.na(dt_all$ens_gene_fruitfly)) |
                                                       (!is.na(dt_all$ens_gene_yeast)) |
                                                       (!is.na(dt_all$ens_gene_mouse)) |
                                                       (!is.na(dt_all$ens_gene_human)) )]
homo_z <- data.table(gene = homo_z$gene, ens_gene = homo_z$ens_gene_zebrafish)

homo_m <- dt_all[(!is.na(dt_all$ens_gene_mouse)) & ( (!is.na(dt_all$ens_gene_fruitfly)) |
                                                       (!is.na(dt_all$ens_gene_zebrafish)) |
                                                       (!is.na(dt_all$ens_gene_yeast)) |
                                                       (!is.na(dt_all$ens_gene_human)) )]
homo_m <- data.table(gene = homo_m$gene, ens_gene = homo_m$ens_gene_mouse)

homo_h <- dt_all[(!is.na(dt_all$ens_gene_human)) & ( (!is.na(dt_all$ens_gene_fruitfly)) |
                                                       (!is.na(dt_all$ens_gene_zebrafish)) |
                                                       (!is.na(dt_all$ens_gene_mouse)) |
                                                       (!is.na(dt_all$ens_gene_yeast)) )]
homo_h <- data.table(gene = homo_h$gene, ens_gene = homo_h$ens_gene_human)


# yeast - sth wrong with IDs...
homo_Y_MF <- enrichGO(gene     = homo_y$ens_gene,
                 OrgDb         = org.Sc.sgd.db,
                 keyType       = 'ENSEMBL', # 'ORF'
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

homo_Y_BP <- enrichGO(gene          = homo_y$ens_gene,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = 'ENSEMBL', # 'ORF'
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# fruit fly
homo_F_MF <- enrichGO(gene          = homo_f$ens_gene,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

homo_F_BP <- enrichGO(gene          = homo_f$ens_gene,
                      OrgDb         = org.Dm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(homo_F_MF, showCategory = 15)
dotplot(homo_F_BP, showCategory = 15)

# zebrafish
homo_Z_MF <- enrichGO(gene          = homo_z$ens_gene,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

homo_Z_BP <- enrichGO(gene          = homo_z$ens_gene,
                      OrgDb         = org.Dr.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(homo_Z_MF, showCategory = 15)
dotplot(homo_Z_BP, showCategory = 15)

# mouse
homo_M_MF <- enrichGO(gene          = homo_m$ens_gene,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

homo_M_BP <- enrichGO(gene          = homo_m$ens_gene,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(homo_M_MF, showCategory = 15)
dotplot(homo_M_BP, showCategory = 15)

# human
homo_H_MF <- enrichGO(gene          = unique(homo_h$ens_gene),
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

homo_H_BP <- enrichGO(gene          = unique(homo_h$ens_gene),
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(homo_H_MF, showCategory = 15)
dotplot(homo_H_BP, showCategory = 15)

#### homologs with stall sites
css <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
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



###################################### load workspace
load("/Volumes/USELESS/STALLING/GO/Rworkspace.RData")

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

wet <- get(load(file = "/Volumes/USELESS/STALLING/peaks_all/wet/human_HeLa_stumpf2013.Rsave"))
wet2 <- get(load(file = "/Volumes/USELESS/STALLING/peaks_all/wet/human_HEK_subtelny2014.Rsave"))
wet3 <- get(load(file = "/Volumes/USELESS/STALLING/peaks_all/wet/human_fibroblasts_stern2012.Rsave")) # has "ENST00000368885"
wet4 <- get(load(file = "/Volumes/USELESS/STALLING/peaks_all/wet/human_stern_NODRUG.Rsave")) # has "ENST00000368885"

wet_tx <- unique(as.character(wet$seqnames))
wet_gene <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx,]$gene_id)
wet_tx2 <- unique(as.character(wet2$seqnames))
wet_gene2 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx2,]$gene_id)
wet_tx3 <- unique(as.character(wet3$seqnames))
wet_gene3 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx3,]$gene_id)
wet_tx4 <- unique(as.character(wet4$seqnames))
wet_gene4 <- unique(txLengths_h[txLengths_h$tx_name %in% wet_tx4,]$gene_id)

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


### define genes without peaks as those with z-score < 3.0 or sth
maxz <- wet[, max(zscore), by = seqnames]
lowz_tx <- as.character(maxz[maxz$V1 < 4]$seqnames)
lowz_gene <- unique(txLengths_h[txLengths_h$tx_name %in% lowz_tx,]$gene_id)
lowz_gene <- lowz_gene[!lowz_gene %in% ss_h]

lowz_H_MF <- enrichGO(gene         = lowz_gene,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

lowz_H_BP <- enrichGO(gene         = lowz_gene,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(lowz_H_MF, showCategory = 30)
dotplot(lowz_H_BP, showCategory = 30)

##########################################################
### export list ss_H_MF and SS_H_BP
### compare to no_wet_H_MF and no_wet_H_BP

ss_H_MF$Description[!ss_H_MF$Description %in% lowz_H_MF$Description]
ss_H_BP$Description[!ss_H_BP$Description %in% lowz_H_BP$Description]

dt_MF_human <- data.table(GO = ss_H_MF$Description, genes = ss_H_MF$geneID)
dt_BP_human <- data.table(GO = ss_H_BP$Description, genes = ss_H_BP$geneID)

write.table(dt_MF_human, "/Volumes/USELESS/STALLING/GO/GO_peaks_human_MF.txt", row.names = F)
write.table(dt_BP_human, "/Volumes/USELESS/STALLING/GO/GO_peaks_human_BP.txt", row.names = F)

### export BED file with human conserved stall sites
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

human_ss <- data.table(tx = css[!is.na(css$h_ss),]$h_tx, ss = css[!is.na(css$h_ss),]$h_ss)
human_ss$end <- human_ss$ss + 2

gr <- makeGRangesFromDataFrame(human_ss,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="tx",
                         start.field="ss",
                         end.field="end",
                         starts.in.df.are.0based=FALSE)

txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
exons_h <- exonsBy(txdb_h, by="tx", use.names=TRUE)

ss_genomic <- mapFromTranscripts(gr, exons_h)

export.bed(ss_genomic, "/Volumes/USELESS/STALLING/GO/GO_peaks_human.bed") ### wrong coord.! nt in cds


