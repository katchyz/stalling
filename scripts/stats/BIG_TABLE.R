### stalling statistics for all genes...
## for human:
## row names: all WET gene/tx names
## columns: #peaks, #cons_peaks, #max_z_score, #GO_enrichment (RNA_metab or not)
library(data.table)
library(ggplot2)

# gene names
homo_path <- "/Volumes/USELESS/STALLING/conservation/homologs"
homo <- get(load(file = c(file.path(homo_path, paste0(org, ".Rsave")))))

go_path <- "/Volumes/USELESS/STALLING/GO/pos_neg_sets"
load(file = c(file.path(go_path, "human_well_expressed_homologs.Rsave"))) # wet_homo, 4860

BIG_TABLE <- data.table(ens_gene = sort(wet_homo))
BIG_TABLE$gene_name <- sapply(BIG_TABLE$ens_gene, function(x){homo[homo$ens_gene == x,]$gene[1]})

# peaks per lib (H1..H4)
h1 <- get(load("/Volumes/USELESS/STALLING/peaks_all/NEW_median/human_h1.Rsave"))
h2 <- get(load("/Volumes/USELESS/STALLING/peaks_all/NEW_median/human_h2.Rsave"))
h3 <- get(load("/Volumes/USELESS/STALLING/peaks_all/NEW_median/human_h3.Rsave"))
h4 <- get(load("/Volumes/USELESS/STALLING/peaks_all/NEW_median/human_h4.Rsave"))

# get unique, count peaks per tx, get highest zscore
all_peaks <- rbind(h1[,1:3], h2[,1:3])
all_peaks <- rbind(all_peaks, h3[,1:3])
all_peaks <- rbind(all_peaks, h4[,1:3])
all_peaks <- unique(all_peaks)

#all_peaks <- all_peaks[all_peaks$seqnames %in% common_human_tx,]
all_peaks$gene <- sapply(as.character(all_peaks$seqnames),
                         function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
load(file = c(file.path(go_path, "human_well_expressed_homologs.Rsave"))) # wet_homo
all_peaks <- all_peaks[all_peaks$gene %in% wet_homo,]

all_peaks$utr5_len <- sapply(as.character(all_peaks$seqnames),
                             function(x){txLengths_h[txLengths_h$tx_name == x,]$utr5_len})
all_peaks$cds_len <- sapply(as.character(all_peaks$seqnames),
                             function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

all_peaks <- all_peaks[all_peaks$start > all_peaks$utr5_len+15]
all_peaks <- all_peaks[all_peaks$start < all_peaks$utr5_len+all_peaks$cds_len-3]

# sum(diff(all_peaks$start) != 3 & diff(all_peaks$start) != 6)
all_peaks <- all_peaks[c(diff(all_peaks$start) != 3, T)] # remove one codon off
all_peaks <- all_peaks[order(seqnames)]

all_peaks$ens_gene <- sapply(as.character(all_peaks$seqnames),
                             function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})

BIG_TABLE$peaks2 <- sapply(BIG_TABLE$ens_gene, function(x){length(all_peaks[all_peaks$gene == x,]$start)})

# css
css_h <- css[!is.na(css$h_tx)]
BIG_TABLE$css <- sapply(BIG_TABLE$gene_name, function(x){length(css_h[css_h$gene == x,]$h_ss)})

# max zscore
# get all 4 libs (exclude 15 and 3nt); then max of max

# wet per lib (H1..H4)
wh1 <- get(load("/Volumes/USELESS/STALLING/peaks_all/wet/human_h1.Rsave"))
wh2 <- get(load("/Volumes/USELESS/STALLING/peaks_all/wet/human_h2.Rsave"))
wh3 <- get(load("/Volumes/USELESS/STALLING/peaks_all/wet/human_h3.Rsave"))
wh4 <- get(load("/Volumes/USELESS/STALLING/peaks_all/wet/human_h4.Rsave"))

#
ltpg <- arrange(txLengths_h, gene_id, desc(tx_len))
ltpg <- ltpg[ltpg$cds_len > 0,]
ltpg <- ltpg[!duplicated(ltpg$gene_id),]
rownames(ltpg) <- ltpg$tx_name

# ltpg
wh1 <- wh1[wh1$seqnames %in% ltpg_h$tx_name,]
wh2 <- wh2[wh2$seqnames %in% ltpg_h$tx_name,]
wh3 <- wh3[wh3$seqnames %in% ltpg_h$tx_name,]
wh4 <- wh4[wh4$seqnames %in% ltpg_h$tx_name,]

# remove first 5 and last 1 for each of seqnames
wh1 <- wh1[, .SD[6:.N-1], by = seqnames]
wh2 <- wh2[, .SD[6:.N-1], by = seqnames]
wh3 <- wh3[, .SD[6:.N-1], by = seqnames]
wh4 <- wh4[, .SD[6:.N-1], by = seqnames]

# remove NA
wh1 <- wh1[!is.na(wh1$zscore)]
wh2 <- wh2[!is.na(wh2$zscore)]
wh3 <- wh3[!is.na(wh3$zscore)]
wh4 <- wh4[!is.na(wh4$zscore)]

mz1 <- wh1[, max(zscore), by = seqnames]
mz2 <- wh2[, max(zscore), by = seqnames]
mz3 <- wh3[, max(zscore), by = seqnames]
mz4 <- wh4[, max(zscore), by = seqnames]

# ltpg
#mz1l <- mz1[mz1$seqnames %in% ltpg_h$tx_name,]
#mz2l <- mz2[mz2$seqnames %in% ltpg_h$tx_name,]
#mz3l <- mz3[mz3$seqnames %in% ltpg_h$tx_name,]
#mz4l <- mz4[mz4$seqnames %in% ltpg_h$tx_name,]

mz <- rbind(mz1, mz2)
mz <- rbind(mz, mz3)
mz <- rbind(mz, mz4)
mz <- mz[order(seqnames)]

mz <- mz[, max(V1), by = seqnames] # min??
mz$gene <- sapply(as.character(mz$seqnames), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
mzg <- mz[, max(V1), by = gene]
#mz1$gene <- sapply(as.character(mz1$seqnames), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
#mz2$gene <- sapply(as.character(mz2$seqnames), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
#mz3$gene <- sapply(as.character(mz3$seqnames), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
#mz4$gene <- sapply(as.character(mz4$seqnames), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})

BIG_TABLE$max_zscore <- as.numeric(sapply(BIG_TABLE$ens_gene, function(x){mz[mz$gene == x,]$V1}))

#BIG_TABLE$maxz_H1 <- sapply(BIG_TABLE$ens_gene, function(x){mz1l[mz1l$gene == x,]$V1})
#BIG_TABLE$maxz_H2 <- sapply(BIG_TABLE$ens_gene, function(x){mz2l[mz2l$gene == x,]$V1})
#BIG_TABLE$maxz_H3 <- sapply(BIG_TABLE$ens_gene, function(x){mz3l[mz3l$gene == x,]$V1})
#BIG_TABLE$maxz_H4 <- sapply(BIG_TABLE$ens_gene, function(x){mz4l[mz4l$gene == x,]$V1})

#BIG_TABLE$maxz_H1 <- as.numeric(BIG_TABLE$maxz_H1)
#BIG_TABLE$maxz_H2 <- as.numeric(BIG_TABLE$maxz_H2)
#BIG_TABLE$maxz_H3 <- as.numeric(BIG_TABLE$maxz_H3)
#BIG_TABLE$maxz_H4 <- as.numeric(BIG_TABLE$maxz_H4)

#BIG_TABLE$maxz <- max(BIG_TABLE$maxz_H1, BIG_TABLE$maxz_H2, BIG_TABLE$maxz_H3, BIG_TABLE$maxz_H4)

#maxz <- c()
#for (i in 1:nrow(BIG_TABLE)) {
#  m <- min(as.numeric(BIG_TABLE[i,5:9])[!is.na(as.numeric(BIG_TABLE[i,5:9]))])
#  maxz <- c(maxz, m)
#}
#BIG_TABLE$maxz <- maxz

# > BIG_TABLE[BIG_TABLE$max_zscore == -Inf]
# ens_gene gene_name peaks css max_zscore
# 1: ENSG00000100814  CCNB1IP1     0   0       -Inf
# 2: ENSG00000104490     NCALD     0   0       -Inf
# 3: ENSG00000110911   SLC11A2     0   0       -Inf
# 4: ENSG00000141255   SPATA22     0   0       -Inf
# 5: ENSG00000171617      ENC1     0   0       -Inf
# 6: ENSG00000181355     OFCC1     0   0       -Inf

## do it only on H2? (so that we get 2064 below zscore 5...)
load(file = "/Volumes/USELESS/STALLING/GO/pos_neg_sets/human_wet_no_peaks.Rsave") # lowz_gene5

## sum(BIG_TABLE$peaks == 0) # 1660 genes with no peaks
## sum(BIG_TABLE[BIG_TABLE$peaks == 0]$ens_gene %in% lowz_gene5) # jak to jest k**wa mozliwe...

# figure out how I got the 2064 lowz_gene5, put it in the table
BIG_TABLE$peaks <- BIG_TABLE$peaks2
BIG_TABLE$peaks2 <- NULL
BIG_TABLE$max_zscore <- NULL
BIG_TABLE$maxz_H1 <- NULL
BIG_TABLE$maxz_H2 <- NULL
BIG_TABLE$maxz_H3 <- NULL
BIG_TABLE$maxz_H4 <- NULL

BIG_TABLE[is.na(BIG_TABLE$max_zscore)]$max_zscore <- BIG_TABLE[is.na(BIG_TABLE$max_zscore)]$maxz

dupa <- sapply(BIG_TABLE[BIG_TABLE$peaks > 0 & BIG_TABLE$max_zscore < 5]$ens_gene,
              function(x){mzg[mzg$gene == x,]$V1})
BIG_TABLE[BIG_TABLE$peaks > 0 & BIG_TABLE$max_zscore < 5]$max_zscore <- as.numeric(dupa)
BIG_TABLE[BIG_TABLE$peaks > 0 & BIG_TABLE$max_zscore < 5]$peaks <- 0
# EMG1 MATR3 

BIG_TABLE[BIG_TABLE$css > 0]$ens_gene # 1730
CSSs$gene <- sapply(as.character(CSSs$tx), function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id})
unique(CSSs$gene) # 1729

BIG_TABLE[BIG_TABLE$css > 0]$ens_gene[!(BIG_TABLE[BIG_TABLE$css > 0]$ens_gene %in% unique(CSSs$gene))]
# "ENSG00000015479" "ENSG00000268439"

BIG_TABLE[BIG_TABLE$ens_gene == "ENSG00000015479"]
BIG_TABLE[BIG_TABLE$ens_gene == "ENSG00000268439"]

BIG_TABLE$maxz <- NULL

### RNA metabolism

BP <- read.table(file = "/Volumes/USELESS/STALLING/GO/WEH_with_peaks_vs_WEH_BP.txt", header = T)
MF <- read.table(file = "/Volumes/USELESS/STALLING/GO/WEH_with_peaks_vs_WEH_MF.txt", header = T)

RNAmetab_genes <- unique(unlist(sapply(BP$genes, function(x){strsplit(as.character(x), split = "/")}))) # 1212
MF_genes <- unique(unlist(sapply(MF$genes, function(x){strsplit(as.character(x), split = "/")})))

BIG_TABLE$rna_metab <- 0
BIG_TABLE[BIG_TABLE$gene_name %in% RNAmetab_genes,]$rna_metab <- 1

save(BIG_TABLE, file = "/Volumes/USELESS/STALLING/BIG_TABLE.Rsave") # 4116 peaks, 738 no peaks
write.csv(BIG_TABLE, file = "/Volumes/USELESS/STALLING/BIG_TABLE.csv")

### PLOT

ggplot(BIG_TABLE[BIG_TABLE$css > 0], aes(x = max_zscore, fill = as.factor(rna_metab))) + geom_histogram(alpha=0.5) + xlim(0,15) + theme_classic() +
  scale_fill_manual(values=c("gray", "#615388"), name = " ", labels = c("all", "RNA_metab"))

ggplot(BIG_TABLE, aes(x = peaks, y = css)) + geom_point()

###

css_MF <- enrichGO(gene          = BIG_TABLE[BIG_TABLE$css > 0]$ens_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   universe      = BIG_TABLE$ens_gene,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

css_BP <- enrichGO(gene          = BIG_TABLE[BIG_TABLE$css > 0]$ens_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   universe      = BIG_TABLE$ens_gene,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(css_MF, showCategory = 25, x = "Count")
dotplot(css_BP, showCategory = 25, x = "Count")

dt_BP <- data.table(GO = css_BP$Description, genes = css_BP$geneID)
unique(unlist(sapply(dt_BP$genes, function(x){strsplit(as.character(x), split = "/")}))) #1212

dt_MF <- data.table(GO = css_MF$Description, genes = css_MF$geneID)
unique(unlist(sapply(dt_MF$genes, function(x){strsplit(as.character(x), split = "/")}))) #766
write.csv(dt_MF, file = "~/Downloads/dt_MF.csv")

### examples
ex <- names(head(sort(table(unlist(sapply(dt_BP$genes,
                                    function(x){strsplit(as.character(x), split = "/")}))), decreasing = T), 10))

css[css$gene %in% ex,]

RP <- RNAmetab_genes[grep("RPS", RNAmetab_genes)]
RP <- c(RP, RNAmetab_genes[grep("RPL", RNAmetab_genes)])
RP_gene <- BIG_TABLE[BIG_TABLE$gene_name %in% RP]$ens_gene
rp_pep <- data.table(tx = css[css$gene %in% RP,]$h_tx, ss_aa = floor(css[css$gene %in% RP,]$h_ss/3)+1)
rp_pep <- rp_pep[complete.cases(rp_pep)]

rp_logo <- list()
for (i in 1:nrow(rp_pep)) {
  tx <- rp_pep[i,]$tx
  ss <- rp_pep[i,]$ss_aa
  seq <- fasta_pep[[tx]][(ss-5):(ss+5)]
  rp_logo <- c(rp_logo, toupper(paste0(seq, collapse = '')))
}

unique(BIG_TABLE$gene_name[grep("RPS", BIG_TABLE$gene_name)]) # 66
unique(BIG_TABLE$gene_name[grep("RPL", BIG_TABLE$gene_name)]) # 99
RPG <- unique(css[grep("RPS", css$gene)]$gene) # 35
RPG <- c(RPG, unique(css[grep("RPL", css$gene)]$gene)) # 53



t <- "ENST00000560055"
ss <- CSSs[CSSs$tx == t]$ss_aa
toupper(paste0(fasta_pep[[t]][(ss-10):(ss+10)], collapse = ''))


