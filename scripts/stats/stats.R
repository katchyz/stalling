### peak stats
library(GenomicFeatures)
library(plyr)
library(data.table)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"

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


# number of peaks and genes with peaks in each library
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
txLengths <- list(txLengths_y, txLengths_f, txLengths_z, txLengths_m, txLengths_h)

for (i in 1:length(organisms)) {
  org <- organisms[i]
  
  ltpg <- arrange(txLengths[i][[1]], gene_id, desc(tx_len))
  ltpg <- ltpg[ltpg$cds_len > 0,]
  ltpg <- ltpg[!duplicated(ltpg$gene_id),]
  rownames(ltpg) <- ltpg$tx_name
  
  libs <- list.files(path = peaks_path, pattern = paste0("^", org))
  
  dupa <- peaks[0,1:3]
  
  for (j in 1:length(libs)) {
    
    print("library:")
    print(libs[j])
    
    ## load file
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- peaks[peaks$seqnames %in% ltpg$tx_name]
    peaks$utr5_len <- sapply(peaks$seqnames, function(x){ltpg[ltpg$tx_name == x,]$utr5_len})
    peaks$cds_len <- sapply(peaks$seqnames, function(x){ltpg[ltpg$tx_name == x,]$cds_len})
    peaks$ss_cds <- peaks$start - peaks$utr5_len
    peaks <- peaks[peaks$ss_cds > 15]
    peaks <- peaks[peaks$ss_cds < peaks$cds_len-3]
    
    #print("# genes")
    #print(length(unique(as.character(peaks$seqnames))))
    #print("# peaks:")
    #print(nrow(peaks))
    
    #print("# peaks in homologs:")
    dupa <- rbind(dupa, peaks[peaks$seqnames %in% homo_human, 1:3])
    
  }
  
  sum(diff(unique(dupa[dupa$seqnames %in% homo_human,])$start) != 3)
}

###################### CONSERVATION ######################
## well-expressed homologs
homo_path <- "/Volumes/USELESS/STALLING/conservation/homologs"
wet_path <- "/Volumes/USELESS/STALLING/peaks_all/wet"
go_path <- "/Volumes/USELESS/STALLING/GO/pos_neg_sets"

css <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
# ss_homo_gene <- unique(ltpg[ltpg$tx_name %in% ss_homo,]$gene)

for (i in 1:length(organisms)) {
  org <- organisms[i]
  
  # ltpg <- arrange(txLengths[i][[1]], gene_id, desc(tx_len))
  # ltpg <- ltpg[ltpg$cds_len > 0,]
  # ltpg <- ltpg[!duplicated(ltpg$gene_id),]
  # rownames(ltpg) <- ltpg$tx_name
  
  # homologs
  homo <- get(load(file = c(file.path(homo_path, paste0(org, ".Rsave")))))
  wet_libs <- list.files(path = wet_path, pattern = paste0("^", org))
  
  wet_all <- get(load(file = c(file.path(wet_path, wet_libs[1]))))
  wet_all <- subset(wet_all,,c(1:7, ncol(wet_all)))
  for (j in 2:length(wet_libs)) {
    ## load file
    wet <- get(load(file = c(file.path(wet_path, wet_libs[j])))) # wet
    wet <- subset(wet,,c(1:7, ncol(wet)))
    wet_all <- rbind(wet_all, wet, fill = T)
  }
  
  wet_tx <- sort(unique(as.character(wet_all$seqnames)))
  wet_gene <- txLengths[i][[1]][txLengths[i][[1]]$tx_name %in% wet_tx,]
  
  print(org)
  print("well expressed genes, total:")
  print(length(unique(wet_gene$gene_id)))
  
  wet_homo <- unique(wet_gene$gene_id)[unique(wet_gene$gene_id) %in% unique(homo$ens_gene)] # all well-expr. homologs
  wet_homo_tx <- unique(wet_gene[wet_gene$gene_id %in% wet_homo,]$tx_name)
  
  print("well expressed homologs:")
  print(length(wet_homo))
  #################
  #save(wet_homo, file = c(file.path(go_path, paste0(org, "_well_expressed_homologs.Rsave"))))
  #################
  
  ss_homo_tx <- unique(as.character(css[!is.na(css[,(i*2+1)]),][,(i*2)]))
  ss_homo_gene <- unique(wet_gene[wet_gene$tx_name %in% ss_homo_tx,]$gene_id)
  
  # ss_homo_gene <- unique(as.character(css[!is.na(css[,(i*2+1)]),]$gene)) ######### gene symbol #########
  
  print("well expressed homologs with peaks:")
  print(length(ss_homo_gene))
  
  wet_all_homo_no <- wet_all[wet_all$seqnames %in% wet_homo_tx] # well expressed homologs
  wet_all_homo_no <- wet_all_homo_no[!wet_all_homo_no$seqnames %in% ss_homo_tx] # well expressed homologs excluding stall sites
  
  maxz <- wet_all_homo_no[, max(zscore), by = seqnames]
  lowz_tx4 <- as.character(maxz[maxz$V1 < 4]$seqnames)
  lowz_gene4 <- unique(txLengths[i][[1]][txLengths[i][[1]]$tx_name %in% lowz_tx4,]$gene_id)
  
  print("well expressed homologs without peaks, zscore < 4:")
  print(length(lowz_gene4))
  
  lowz_tx5 <- as.character(maxz[maxz$V1 < 5]$seqnames)
  lowz_gene5 <- unique(txLengths[i][[1]][txLengths[i][[1]]$tx_name %in% lowz_tx5,]$gene_id)
  
  save(lowz_gene5, file = c(file.path(go_path, paste0(org, "_wet_no_peaks.Rsave"))))
  
  print("well expressed homologs without peaks, zscore < 5:")
  print(length(lowz_gene5))
  
}


# for human:
save(ss_homo_gene, file = c(file.path(go_path, paste0(org, "_homologs_with_peaks_gt5zscore.Rsave"))))
save(lowz_gene4, file = c(file.path(go_path, paste0(org, "_homologs_without_peaks_lt4zscore.Rsave"))))
save(lowz_gene5, file = c(file.path(go_path, paste0(org, "_homologs_without_peaks_lt4zscore.Rsave"))))
save(wet_homo, file = c(file.path(go_path, paste0(org, "_well_expressed_homologs.Rsave"))))


##########
# check GO
library(clusterProfiler)
library(org.Hs.eg.db)

ss_H_MF <- enrichGO(gene       = ss_homo_gene,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

ss_H_BP <- enrichGO(gene       = ss_homo_gene,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

ss_H_CC <- enrichGO(gene       = ss_homo_gene,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_H_MF4 <- enrichGO(gene          = lowz_gene4,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_H_BP4 <- enrichGO(gene          = lowz_gene4,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

no_H_CC4 <- enrichGO(gene          = lowz_gene4,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

no_H_MF5 <- enrichGO(gene          = lowz_gene5,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

no_H_BP5 <- enrichGO(gene          = lowz_gene5,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

no_H_CC5 <- enrichGO(gene          = lowz_gene5,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(ss_H_MF, showCategory = 30, x = "Count")
dotplot(ss_H_BP, showCategory = 30, x = "Count")
# dotplot(ss_H_CC, showCategory = 30, x = "Count")
dotplot(no_H_MF4, showCategory = 30, x = "Count")
dotplot(no_H_BP4, showCategory = 30, x = "Count")
# dotplot(no_H_CC4, showCategory = 30, x = "Count")
dotplot(no_H_MF5, showCategory = 30, x = "Count")
dotplot(no_H_BP5, showCategory = 30, x = "Count")
# dotplot(no_H_CC5, showCategory = 30, x = "Count")

#######

bckg_H_MF <- enrichGO(gene          = ss_homo_gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = wet_homo,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_H_BP <- enrichGO(gene          = ss_homo_gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = wet_homo,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_H_CC <- enrichGO(gene          = ss_homo_gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      universe      = wet_homo,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_no_H_MF4 <- enrichGO(gene          = lowz_gene4,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      universe      = wet_homo,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_no_H_BP4 <- enrichGO(gene          = lowz_gene4,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      universe      = wet_homo,
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

bckg_no_H_CC4 <- enrichGO(gene          = lowz_gene4,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC",
                          pAdjustMethod = "BH",
                          universe      = wet_homo,
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

bckg_no_H_MF5 <- enrichGO(gene          = lowz_gene5,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         universe      = wet_homo,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)

bckg_no_H_BP5 <- enrichGO(gene          = lowz_gene5,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         universe      = wet_homo,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)

bckg_no_H_CC5 <- enrichGO(gene          = lowz_gene5,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC",
                          pAdjustMethod = "BH",
                          universe      = wet_homo,
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

dotplot(bckg_H_MF, showCategory = 45, x = "Count")
dotplot(bckg_H_BP, showCategory = 45, x = "Count")
dotplot(bckg_H_CC, showCategory = 45, x = "Count")
dotplot(bckg_no_H_MF4, showCategory = 30, x = "Count")
dotplot(bckg_no_H_BP4, showCategory = 30, x = "Count")
dotplot(bckg_no_H_CC4, showCategory = 30, x = "Count")
dotplot(bckg_no_H_MF5, showCategory = 30, x = "Count")
dotplot(bckg_no_H_BP5, showCategory = 30, x = "Count")
dotplot(bckg_no_H_CC5, showCategory = 30, x = "Count")
dt_BP_human <- data.table(GO = bckg_H_BP$Description, genes = bckg_H_BP$geneID)
write.table(dt_BP_human, "/Volumes/USELESS/STALLING/GO/WEH_with_peaks_vs_WEH_BP.txt", row.names = F)
dt_MF_human <- data.table(GO = bckg_H_MF$Description, genes = bckg_H_MF$geneID)
write.table(dt_MF_human, "/Volumes/USELESS/STALLING/GO/WEH_with_peaks_vs_WEH_MF.txt", row.names = F)
dt_CC_human <- data.table(GO = bckg_H_CC$Description, genes = bckg_H_CC$geneID)
write.table(dt_CC_human, "/Volumes/USELESS/STALLING/GO/WEH_with_peaks_vs_WEH_CC.txt", row.names = F)


# dt_BP_human <- data.table(GO = ss_H_BP$Description, genes = ss_H_BP$geneID)
# write.table(dt_BP_human, "/Volumes/USELESS/STALLING/GO/UPDATED_GO_peaks_human_BP.txt", row.names = F)
# dt_MF_human <- data.table(GO = ss_H_MF$Description, genes = ss_H_MF$geneID)
# write.table(dt_MF_human, "/Volumes/USELESS/STALLING/GO/UPDATED_GO_peaks_human_MF.txt", row.names = F)

# fileConn<-file("/Volumes/USELESS/STALLING/GO/homologs_with_stall_sites.txt")
# writeLines(ss_homo_gene, fileConn)
# close(fileConn)

# ENST00000216037
# library(ggplot2)
# xbp1 <- wet_all[wet_all$seqnames == "ENST00000216037"][1:259]
# some_gene <- wet_all[wet_all$seqnames == "ENST00000009180"]
# ggplot(some_gene, aes(x = zscore)) + geom_density() #+ theme_void()


# ENST00000005259 ENST00000007516 ENST00000009180 ENST00000009589 ENST00000010404 ENST00000011473 ENST00000011691 ENST00000011898


ggplot(bckg_H_BP, aes(x = Count, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")





