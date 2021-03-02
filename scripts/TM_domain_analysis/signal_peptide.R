## signal peptide

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]


signal_peptide_path <- "../../DATA/domains"
save_plot_path <- "../../DATA/domains"
# TM1 #############
sp_files <- c("tm1_yeast_gene_name.txt", "tm1_fruitfly_flybase.txt", "tm1_zebrafish_gene.txt", "tm1_mouse_gene.txt", "tm1_human_gene.txt")
# TM2 #############
sp_files <- c("tm2_yeast_gene_name.txt", "tm2_fruitfly_flybase.txt", "tm2_zebrafish_gene.txt", "tm2_mouse_gene.txt", "tm2_human_gene.txt")

######################
### CSSs
org <- "human"
css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
#
sp <- as.character(read.table(c(file.path(signal_peptide_path, sp_files[5])))$V1)
#
sp_tx <- txLengths_h[txLengths_h$gene_id %in% sp,]$tx_name
consensus <- css[!is.na(css$h_ss),]
consensus <- consensus[consensus$h_tx %in% sp_tx,]
consensus$aa <- floor(consensus$h_ss/3) + 1
print(nrow(consensus))
## plot from start
ggplot(consensus, aes(x=aa)) + geom_histogram(binwidth = 1) + xlim(0,100)
# ggsave(file = c(file.path(save_plot_path, paste0(org, ".png"))))


