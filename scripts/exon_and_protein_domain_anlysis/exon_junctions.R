### distance to exon junctions

# metaplot around junctions
# optionally: each exon normalized to 1, metaplot

### HUMAN
library(GenomicFeatures)
library(data.table)
library(ggplot2)

# load genome, get exons
gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]
cds_h <- cdsBy(txdb_h, by="tx", use.names=TRUE)

# get stall site positions (on CDS)
css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
ss_h <- data.table(tx = as.character(css[!is.na(css$h_ss),]$h_tx), ss = css[!is.na(css$h_ss),]$h_ss)

# subset tx with peaks
cds <- cds_h[names(cds_h) %in% ss_h$tx] # 1729 tx

gr <- makeGRangesFromDataFrame(ss_h,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="tx",
                               start.field="ss",
                               end.field="ss",
                               starts.in.df.are.0based=FALSE)

ss_genomic <- mapFromTranscripts(gr, cds_h)
ss_genomic$seqnames <- as.character(seqnames(gr))
# export.bed(ss_genomic, "../../DATA/GO/GO_peaks_human.bed")

# cross-check exons with stall site positions

# what exon each stall site belong to
fover2 <- findOverlaps(cds, ss_genomic)
start <- rep(NA, length(ss_genomic))
start[to(fover2)] <- start(cds)[from(fover2)]
end <- rep(NA, length(ss_genomic))
end[to(fover2)] <- end(cds)[from(fover2)]
elementMetadata(ss_genomic)$exons <- IRanges(start=start, end=end)

# calculate distance to junction before and junction after
plus <- as.character(strand(ss_genomic)) == "+"
minus <- as.character(strand(ss_genomic)) == "-"
junction_before <- rep(NA, length(ss_genomic))
junction_before[plus] <- start(ss_genomic)[plus] - start(ss_genomic$exons)[plus]
junction_before[minus] <- end(ss_genomic$exons)[minus] - start(ss_genomic)[minus]
junction_after <- rep(NA, length(ss_genomic))
junction_after[plus] <- start(ss_genomic)[plus] - end(ss_genomic$exons)[plus]
junction_after[minus] <- start(ss_genomic$exons)[minus] - start(ss_genomic)[minus]

elementMetadata(ss_genomic)$junction_before <- junction_before
elementMetadata(ss_genomic)$junction_after <- junction_after

elementMetadata(ss_genomic)$exon_width <- width(ss_genomic$exons)

dt_ss <- as.data.table(ss_genomic)
dt_ss$position_on_exon <- dt_ss$junction_before / dt_ss$exon_width

ggplot(dt_ss, aes(x = position_on_exon)) + geom_histogram()
# ggsave("position_on_exons.png")

dt_junc <- data.table(type = c(rep("junction_before", nrow(dt_ss)), rep("junction_after", nrow(dt_ss))),
                      distance = c(dt_ss$junction_before, dt_ss$junction_after))

ggplot(dt_junc, aes(x = distance)) + geom_histogram(binwidth = 1) + xlim(-50,51)
# ggsave("distance_to_junction_bin1.png")

sort(table(dt_junc$distance), decreasing = T) # most common: 14, 20, 17, 23, 29, 26
# frames:
# 0    1    2 
# 496  640 1290

#### which exon (is it last one?)
# how many stall sites at given exon (check exon rank)
cds <- unlist(cds)
fover <- findOverlaps(ss_genomic, cds)
cds$ss <- rep(NA, length(cds))
cds$ss[as.numeric(names(table(to(fover))))] <- table(to(fover))

cds$tx_name <- names(cds)
dt_cds <- as.data.table(cds)
dt_cds <- dt_cds[order(tx_name, start)]

dt_cds$exon_from_last <- unlist(sapply(table(dt_cds$tx_name), function(x){c(-x:-1)}))
ggplot(dt_cds, aes(x = exon_from_last)) + geom_histogram() + xlim(-10,0)

########################################################################################################
####### get subsets #######
# load GO terms (BP)
go <- read.table("../../DATA/GO/GO_peaks_human_BP.txt", header = T)
splicing_genes <- as.character(go$genes[grep("splicing", go$GO)])
splicing_genes <- unique(unlist(sapply(splicing_genes, function(x){strsplit(x, split = "/")})))
splicing_tx <- as.character(css[css$gene %in% splicing_genes & !is.na(css$h_ss),]$h_tx)

splicing_ss <- dt_ss[dt_ss$seqnames.1 %in% splicing_tx]
ggplot(splicing_ss, aes(x = position_on_exon)) + geom_histogram()

splicing_junc <- data.table(type = c(rep("junction_before", nrow(splicing_ss)), rep("junction_after", nrow(splicing_ss))),
                      distance = c(splicing_ss$junction_before, splicing_ss$junction_after))
ggplot(splicing_junc, aes(x = distance)) + geom_histogram(binwidth = 1) + xlim(-50,51)
# ggsave("distance_to_junction_SPLICING_GO_bin1.png")

# splicing_genes
splicing14_tx <- splicing_ss[splicing_ss$junction_before == 14]$seqnames.1
splicing14_gene <- unique(css[css$h_tx %in% splicing14_tx,]$gene)

# 41 stall sites in 17 genes

# GO on the ones stalled at distance of 14 from exon junction
all14_tx <- dt_ss[dt_ss$junction_before == 14]$seqnames.1
all14_gene <- unique(sapply(all14_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$gene_id}))

library(clusterProfiler)
library(org.Hs.eg.db)

h14_CC <- enrichGO(gene          = all14_gene,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

h14_MF <- enrichGO(gene          = all14_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

h14_BP <- enrichGO(gene          = all14_gene,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(h14_CC)
dotplot(h14_MF)
dotplot(h14_BP)

########
# load("R_exon_junctions.RData")

# check if those 14nt downstream are in 1st exon as well
# in other words: check if the stall site position is 14 as well
seq14 <- ss_genomic[ss_genomic$junction_before == 14]$seqnames
sum(css[css$h_tx %in% seq14,]$h_ss == 15)


# ignore first CDS exon
no1ex_plus <- cds[strand(cds) == "+"]
no1ex_plus <- split(no1ex_plus, names(no1ex_plus))
no1ex_plus <- no1ex_plus[sapply(no1ex_plus, function(x){length(x)}) > 1]
no1ex_minus <- cds[strand(cds) == "-"]
no1ex_minus <- split(no1ex_minus, names(no1ex_minus))
no1ex_minus <- no1ex_minus[sapply(no1ex_minus, function(x){length(x)}) > 1]

n1p <- GRangesList(sapply(no1ex_plus, function(x){x[2:length(x)]}))
n1m <- GRangesList(sapply(no1ex_minus, function(x){x[2:length(x)]}))

no1ex <- c(n1p, n1m)

# what exon each stall site belong to
no1ex <- unlist(no1ex)
fover3 <- findOverlaps(no1ex, ss_genomic)
start <- rep(NA, length(ss_genomic))
start[to(fover3)] <- start(no1ex)[from(fover3)]
end <- rep(NA, length(ss_genomic))
end[to(fover3)] <- end(no1ex)[from(fover3)]
ss_genomic <- ss_genomic[!is.na(start)]
start <- start[!is.na(start)]
end <- end[!is.na(end)]
elementMetadata(ss_genomic)$exons <- IRanges(start=start, end=end)

# calculate distance to junction before and junction after
plus <- as.character(strand(ss_genomic)) == "+"
minus <- as.character(strand(ss_genomic)) == "-"
junction_before <- rep(NA, length(ss_genomic))
junction_before[plus] <- start(ss_genomic)[plus] - start(ss_genomic$exons)[plus]
junction_before[minus] <- end(ss_genomic$exons)[minus] - start(ss_genomic)[minus]
junction_after <- rep(NA, length(ss_genomic))
junction_after[plus] <- start(ss_genomic)[plus] - end(ss_genomic$exons)[plus]
junction_after[minus] <- start(ss_genomic$exons)[minus] - start(ss_genomic)[minus]

elementMetadata(ss_genomic)$junction_before <- junction_before
elementMetadata(ss_genomic)$junction_after <- junction_after

elementMetadata(ss_genomic)$exon_width <- width(ss_genomic$exons)

dt_ss <- as.data.table(ss_genomic)
dt_ss$position_on_exon <- dt_ss$junction_before / dt_ss$exon_width

ggplot(dt_ss, aes(x = position_on_exon)) + geom_histogram()
# ggsave("position_on_exons_NO1EX.png")

dt_junc <- data.table(type = c(rep("junction_before", nrow(dt_ss)), rep("junction_after", nrow(dt_ss))),
                      distance = c(dt_ss$junction_before, dt_ss$junction_after))

ggplot(dt_junc, aes(x = distance)) + geom_histogram(binwidth = 1) + xlim(-50,51)
# ggsave("distance_to_junction_bin1_NO1EX.png")


#################################################################################
#################################################################################
### redo, remove first 15 CODONS (stall site > 45 on CDS)

# get stall site positions (on CDS)
ss_h <- ss_h[ss_h$ss > 45] # 2030 genes

# subset tx with peaks
cds <- cds_h[names(cds_h) %in% ss_h$tx] # 1729 tx

gr <- makeGRangesFromDataFrame(ss_h,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="tx",
                               start.field="ss",
                               end.field="ss",
                               starts.in.df.are.0based=FALSE)

ss_genomic <- mapFromTranscripts(gr, cds_h)
ss_genomic$seqnames <- as.character(seqnames(gr))
# export.bed(ss_genomic, "../../DATA/GO/GO_peaks_human.bed")

# cross-check exons with stall site positions

# what exon each stall site belong to
cds <- unlist(cds)
fover2 <- findOverlaps(cds, ss_genomic)
start <- rep(NA, length(ss_genomic))
start[to(fover2)] <- start(cds)[from(fover2)]
end <- rep(NA, length(ss_genomic))
end[to(fover2)] <- end(cds)[from(fover2)]
elementMetadata(ss_genomic)$exons <- IRanges(start=start, end=end)

# calculate distance to junction before and junction after
plus <- as.character(strand(ss_genomic)) == "+"
minus <- as.character(strand(ss_genomic)) == "-"
junction_before <- rep(NA, length(ss_genomic))
junction_before[plus] <- start(ss_genomic)[plus] - start(ss_genomic$exons)[plus]
junction_before[minus] <- end(ss_genomic$exons)[minus] - start(ss_genomic)[minus]
junction_after <- rep(NA, length(ss_genomic))
junction_after[plus] <- start(ss_genomic)[plus] - end(ss_genomic$exons)[plus]
junction_after[minus] <- start(ss_genomic$exons)[minus] - start(ss_genomic)[minus]

elementMetadata(ss_genomic)$junction_before <- junction_before
elementMetadata(ss_genomic)$junction_after <- junction_after

elementMetadata(ss_genomic)$exon_width <- width(ss_genomic$exons)

dt_ss <- as.data.table(ss_genomic)
dt_ss$position_on_exon <- dt_ss$junction_before / dt_ss$exon_width

ggplot(dt_ss, aes(x = position_on_exon)) + geom_histogram() + theme_classic()
# ggsave("position_on_exons_start45.png")
# save(dt_ss, file = "../../DATA/domains/CSSs_position_on_exon.Rsave")

dt_junc <- data.table(type = c(rep("junction_before", nrow(dt_ss)), rep("junction_after", nrow(dt_ss))),
                      distance = c(dt_ss$junction_before, dt_ss$junction_after))


ggplot(dt_junc, aes(x = distance)) + geom_histogram(binwidth = 1) + xlim(-50,51) + theme_classic()
# ggsave("distance_to_junction_bin1_start45.png")
# change direction
dt_junc$distance <- -dt_junc$distance
ggplot(dt_junc, aes(x = distance)) + geom_histogram(binwidth = 1) + xlim(-50,51) + theme_classic()
# ggsave("distance_to_junction_bin1_start45_reversed.png")

sort(table(dt_junc$distance), decreasing = T) # most common: -15, 0, 8, -9, -18, -1



