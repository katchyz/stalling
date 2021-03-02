library(data.table)
library(GenomicFeatures)
library(seqinr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Sc.sgd.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)

css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
setDT(css)

### conserved in 3 organisms (incl. human)

cons3 <- css[as.numeric(!is.na(css$y_tx)) +
               as.numeric(!is.na(css$f_tx)) +
               as.numeric(!is.na(css$z_tx)) +
               as.numeric(!is.na(css$m_tx)) +
               as.numeric(!is.na(css$h_tx)) > 2]

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]

cons3_h <- cons3[!is.na(cons3$h_tx)]
#cons3_h$utr5_len <- sapply(as.character(cons3_h$h_tx), function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})

c3h <- data.table(tx = cons3_h$h_tx, ss_aa = floor(cons3_h$h_ss/3)+1 )
c3h <- c3h[order(c3h$tx),]

save(c3h, file = "../../DATA/sequence/cons3_human.Rsave")

##
load("../../DATA/sequence/fasta_pep.Rsave")

peps <- c()
for (t in as.character(c3h$tx)) {
  ss <- c3h[tx == t]$ss_aa
  pep <- fasta_pep[[t]][(ss-5):(ss+5)]
  peps <- c(peps, paste0(toupper(pep), collapse = ''))
}

#pep <- mapply(function(x,y){fasta_pep[[x]][(y-5):(y+5)]}, as.character(c3h$tx), c3h$ss_aa)
#sapply(pep, function(x){paste0(toupper(x), collapse = '')})

c3_tx <- unique(as.character(cons3$h_tx[!is.na(cons3$h_tx)]))
c3_genes <- sapply(c3_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})

mouse <- data.table(tx = mouse$m_tx, ss_aa = floor(mouse$m_ss/3)+1)

go_path <- "../../DATA/GO"
load(file = c(file.path(go_path, "human_well_expressed_homologs.Rsave"))) # wet_homo

cons3_MF <- enrichGO(gene          = c3_genes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wet_homo,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

cons3_BP <- enrichGO(gene          = c3_genes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wet_homo,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(cons3_MF, showCategory = 35, x = "Count")
dotplot(cons3_BP, showCategory = 35, x = "Count")

dt_BP_cons3 <- data.table(GO = cons3_BP$Description, genes = cons3_BP$geneID)
#write.table(dt_BP_cons3, "../../DATA/GO/cons3_BP.txt", row.names = F)
#dt_MF_cons3 <- data.table(GO = cons3_MF$Description, genes = cons3_MF$geneID)
#write.table(dt_MF_cons3, "../../DATA/GO/cons3_MF.txt", row.names = F)


#### repeat for zebrafish, fruit fly, yeast

## zebrafish
mouse <- css[!is.na(css$m_tx)]
zebra <- css[!is.na(css$z_tx)]
fruit <- css[!is.na(css$f_tx)]
yeast <- css[!is.na(css$y_tx)]

mouse <- data.table(tx = mouse$m_tx, ss_aa = floor(mouse$m_ss/3)+1)
zebra <- data.table(tx = zebra$z_tx, ss_aa = floor(zebra$z_ss/3)+1)
fruit <- data.table(tx = fruit$f_tx, ss_aa = floor(fruit$f_ss/3)+1)
yeast <- data.table(tx = yeast$y_tx, ss_aa = floor(yeast$y_ss/3)+1)

#write.csv(mouse, file = "../../DATA/sequence/mouse.txt", row.names = F, quote = F)
#write.csv(zebra, file = "../../DATA/sequence/zebra.txt", row.names = F, quote = F)
#write.csv(fruit, file = "../../DATA/sequence/fruit.txt", row.names = F, quote = F)
#write.csv(yeast, file = "../../DATA/sequence/yeast.txt", row.names = F, quote = F)
#write.csv(c3h, file = "../../DATA/sequence/cons3_human.txt", row.names = F, quote = F)

# fasta for AA analysis
fasta_yeast <- read.fasta("../../DATA/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz")
fasta_fruit <- read.fasta("../../DATA/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.pep.all.fa.gz")
fasta_zebra <- read.fasta("../../DATA/fasta/danio_rerio/GRCz10/pep/Danio_rerio.GRCz10.pep.all.fa.gz")
fasta_mouse <- read.fasta("../../DATA/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz")

names(fasta_fruit) <- sapply(getAnnot(fasta_fruit), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})
names(fasta_zebra) <- sapply(getAnnot(fasta_zebra), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})
names(fasta_zebra) <- sapply(names(fasta_zebra), function(x){substr(x,1,18)})
names(fasta_mouse) <- sapply(getAnnot(fasta_mouse), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})

fasta_yeast <- fasta_yeast[names(fasta_yeast) %in% yeast$y_tx]
fasta_fruit <- fasta_fruit[names(fasta_fruit) %in% fruit$f_tx]
fasta_zebra <- fasta_zebra[names(fasta_zebra) %in% zebra$z_tx]
fasta_mouse <- fasta_mouse[names(fasta_mouse) %in% mouse$tx]
fasta_cons3 <- fasta_pep[names(fasta_pep) %in% c3h$tx]

#write.fasta(sequences = fasta_yeast, names = names(fasta_yeast),
#            file.out = "../../DATA/sequence/fasta_yeast.fa")
#write.fasta(sequences = fasta_fruit, names = names(fasta_fruit),
#            file.out = "../../DATA/sequence/fasta_fruit.fa")
#write.fasta(sequences = fasta_zebra, names = names(fasta_zebra),
#            file.out = "../../DATA/sequence/fasta_zebra.fa")
#write.fasta(sequences = fasta_mouse, names = names(fasta_mouse),
#            file.out = "../../DATA/sequence/fasta_mouse.fa")
#write.fasta(sequences = fasta_cons3, names = names(fasta_cons3),
#            file.out = "../../DATA/sequence/fasta_cons3_human.fa")

# txdb and wet_homo for GO analysis
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

yeast$gene <- sapply(as.character(yeast$tx), function(x){txLengths_y[txLengths_y$tx_name == x,]$gene_id})
fruit$gene <- sapply(as.character(fruit$tx), function(x){txLengths_f[txLengths_f$tx_name == x,]$gene_id})
zebra$gene <- sapply(as.character(zebra$tx), function(x){txLengths_z[txLengths_z$tx_name == x,]$gene_id})
mouse$gene <- sapply(as.character(mouse$tx), function(x){txLengths_m[txLengths_m$tx_name == x,]$gene_id})

wethomo_yeast <- get(load(file = c(file.path(go_path, "yeast_well_expressed_homologs.Rsave"))))
wethomo_fruit <- get(load(file = c(file.path(go_path, "fruitfly_well_expressed_homologs.Rsave"))))
wethomo_zebra <- get(load(file = c(file.path(go_path, "zebrafish_well_expressed_homologs.Rsave"))))
wethomo_mouse <- get(load(file = c(file.path(go_path, "mouse_well_expressed_homologs.Rsave"))))

## GO

# zebra
zebra_MF <- enrichGO(gene          = zebra$gene,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_zebra,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

zebra_BP <- enrichGO(gene          = zebra$gene,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_zebra,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(zebra_MF, showCategory = 25, x = "Count")
dotplot(zebra_BP, showCategory = 15, x = "Count")

dt_BP_zebra <- data.table(GO = zebra_BP$Description, genes = zebra_BP$geneID)
# write.table(dt_BP_zebra, "../../DATA/GO/zebra_BP.txt", row.names = F)
dt_MF_zebra <- data.table(GO = zebra_MF$Description, genes = zebra_MF$geneID)
# write.table(dt_MF_zebra, "../../DATA/GO/zebra_MF.txt", row.names = F)

# fruit
fruit_MF <- enrichGO(gene          = fruit$gene,
                     OrgDb         = org.Dm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_fruit,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

fruit_BP <- enrichGO(gene          = fruit$gene,
                     OrgDb         = org.Dm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_fruit,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(fruit_MF, showCategory = 15, x = "Count")
dotplot(fruit_BP, showCategory = 15, x = "Count")

dt_BP_fruit <- data.table(GO = fruit_BP$Description, genes = fruit_BP$geneID)
# write.table(dt_BP_fruit, ".././DATA/GO/fruit_BP.txt", row.names = F)
dt_MF_fruit <- data.table(GO = fruit_MF$Description, genes = fruit_MF$geneID)
# write.table(dt_MF_fruit, "../../DATA/GO/fruit_MF.txt", row.names = F)

# yeast
#yeast_entrez <- bitr(yeast$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Sc.sgd.db)
#yeast_wethomo_entrez <- bitr(wethomo_yeast, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Sc.sgd.db)

yeast_MF <- enrichGO(gene          = yeast$gene,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_yeast,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

yeast_BP <- enrichGO(gene          = yeast$gene,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_yeast,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

dotplot(yeast_MF, showCategory = 20, x = "Count")
dotplot(yeast_BP, showCategory = 15, x = "Count")

dt_BP_yeast <- data.table(GO = yeast_BP$Description, genes = yeast_BP$geneID)
# write.table(dt_BP_yeast, "../../DATA/GO/yeast_BP.txt", row.names = F)
dt_MF_yeast <- data.table(GO = yeast_MF$Description, genes = yeast_MF$geneID)
# write.table(dt_MF_yeast, "../../DATA/GO/yeast_MF.txt", row.names = F)


# mouse
mouse_MF <- enrichGO(gene          = mouse$gene,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_mouse,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

mouse_BP <- enrichGO(gene          = mouse$gene,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_mouse,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(mouse_MF, showCategory = 25, x = "Count")
dotplot(mouse_BP, showCategory = 25, x = "Count")

dt_BP_mouse <- data.table(GO = mouse_BP$Description, genes = mouse_BP$geneID)
# write.table(dt_BP_mouse, "../../DATA/GO/mouse_BP.txt", row.names = F)
dt_MF_mouse <- data.table(GO = mouse_MF$Description, genes = mouse_MF$geneID)
# write.table(dt_MF_mouse, "../../DATA/GO/mouse_MF.txt", row.names = F)

### control!
go_path <- "../../DATA/GO"
ctr_yeast <- get(load(file = c(file.path(go_path, "yeast_wet_no_peaks.Rsave")))) #14
ctr_fruit <- get(load(file = c(file.path(go_path, "fruitfly_wet_no_peaks.Rsave")))) #68
ctr_zebra <- get(load(file = c(file.path(go_path, "zebrafish_wet_no_peaks.Rsave")))) #362
ctr_mouse <- get(load(file = c(file.path(go_path, "mouse_wet_no_peaks.Rsave")))) #1373

# mouse
mouse_ctr_MF <- enrichGO(gene      = ctr_mouse,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_mouse,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

mouse_ctr_BP <- enrichGO(gene      = ctr_mouse,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_mouse,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(mouse_ctr_MF, showCategory = 25, x = "Count")
dotplot(mouse_ctr_BP, showCategory = 25, x = "Count")

# zebra
zebra_ctr_MF <- enrichGO(gene      = ctr_zebra,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_zebra,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

zebra_ctr_BP <- enrichGO(gene          = ctr_zebra,
                     OrgDb         = org.Dr.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_zebra,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(zebra_ctr_MF, showCategory = 25, x = "Count")
dotplot(zebra_ctr_BP, showCategory = 15, x = "Count")

# fruit
fruit_ctr_MF <- enrichGO(gene      = ctr_fruit,
                     OrgDb         = org.Dm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_fruit,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

fruit_ctr_BP <- enrichGO(gene      = ctr_fruit,
                     OrgDb         = org.Dm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_fruit,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(fruit_ctr_MF, showCategory = 15, x = "Count")
dotplot(fruit_ctr_BP, showCategory = 15, x = "Count")

yeast_ctr_MF <- enrichGO(gene      = ctr_yeast,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = 'ENSEMBL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     universe      = wethomo_yeast,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

yeast_ctr_BP <- enrichGO(gene      = ctr_yeast,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     universe      = wethomo_yeast,
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

dotplot(yeast_ctr_MF, showCategory = 20, x = "Count")
dotplot(yeast_ctr_BP, showCategory = 15, x = "Count")


