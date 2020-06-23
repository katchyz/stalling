data1 <- read.csv('/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median_all_peaks/NEW_median_all_peaks.csv', header = TRUE)

#MF <- read.table("/Volumes/USELESS/STALLING/GO/GO_peaks_human_MF.txt")
#BP <- read.table("/Volumes/USELESS/STALLING/GO/GO_peaks_human_BP.txt")

BP <- read.table("/Users/kasia/Desktop/dupa/GO.csv", header = T)
BP$X <- NULL

go <- sapply(as.character(BP$genes), function(x){strsplit(x, split = "/")})
names(go) <- BP$GO

data1$gene <- sapply(as.character(data1$gene), function(x){strsplit(x, split = "_")[[1]][1]})

hm <- data1[data1$human == 1 & data1$mouse == 1 & data1$zebrafish == 0 & data1$fruitfly == 0 & data1$yeast == 0,]$gene
nohuman <- unique(data1[data1$human == 0 & data1$mouse == 1,]$gene)
nohumannomouse <- unique(data1[data1$human == 0 & data1$mouse == 0,]$gene)

css_nohuman <- css[!is.na(css$m_ss) & is.na(css$h_ss),]
css_nohumannomouse <- css[is.na(css$m_ss) & is.na(css$h_ss),]
nohuman_ensembl <- unique(txLengths[txLengths$tx_name %in% css_nohuman$m_tx,]$gene_id)
nohumannomouse_ensembl <- unique(txLengths_z[txLengths_z$tx_name %in% css_nohumannomouse$z_tx,]$gene_id)


nohuman_BP <- enrichGO(gene          = nohuman_ensembl,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(nohuman_BP, showCategory = 30, x = "Count")

nohuman_MF <- enrichGO(gene          = nohuman_ensembl,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

dotplot(nohuman_MF, showCategory = 30, x = "Count")

nohumannomouse_BP <- enrichGO(gene          = nohumannomouse_ensembl,
                       OrgDb         = org.Dr.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

dotplot(nohumannomouse_BP, showCategory = 30, x = "Count")

nohumannomouse_MF <- enrichGO(gene          = nohumannomouse_ensembl,
                              OrgDb         = org.Dr.eg.db,
                              keyType       = 'ENSEMBL',
                              ont           = "MF",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

dotplot(nohumannomouse_MF, showCategory = 30, x = "Count")




