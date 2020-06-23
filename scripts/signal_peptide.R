## signal peptide

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
signal_peptide_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/signal_peptide_proteins"
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/TM"
sp_files <- c("yeast_sp_gene_name.txt", "fruitfly_sp_flybase.txt", "zebrafish_sp_gene.txt", "mouse_sp_gene.txt", "human_sp_gene.txt")
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

for (i in 1:length(organisms)) {
  org <- organisms[i]
  sp <- as.character(read.table(c(file.path(signal_peptide_path, sp_files[i])))$V1)
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave"))))
  #
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  txLengths <- txLengths[order(txLengths$tx_name),]
  #
  sp_tx <- txLengths[txLengths$gene_id %in% sp,]$tx_name
  print(org)
  consensus <- consensus[consensus$seqnames %in% sp_tx,]
  print(nrow(consensus))
  ## plot from start
  ggplot(consensus, aes(x=ss)) + geom_histogram(bins = nrow(consensus))
  ggsave(file = c(file.path(save_plot_path, paste0(org, ".png"))))
  ggplot(consensus, aes(x=ss)) + geom_histogram() + scale_x_log10()
  ggsave(file = c(file.path(save_plot_path, paste0(org, "_log10.png"))))
}


############# TM1 #############
# consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
signal_peptide_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/tm1"
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/tm1_plots"
sp_files <- c("tm1_yeast_gene_name.txt", "tm1_fruitfly_flybase.txt", "tm1_zebrafish_gene.txt", "tm1_mouse_gene.txt", "tm1_human_gene.txt")
# organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

############# TM2 #############

# consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
signal_peptide_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/tm2"
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/tm2_plots"
sp_files <- c("tm2_yeast_gene_name.txt", "tm2_fruitfly_flybase.txt", "tm2_zebrafish_gene.txt", "tm2_mouse_gene.txt", "tm2_human_gene.txt")
# organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")



######################
### CSSs
org <- "human"
sp <- as.character(read.table(c(file.path(signal_peptide_path, sp_files[5])))$V1)
#
sp_tx <- txLengths_h[txLengths_h$gene_id %in% sp,]$tx_name
consensus <- css[!is.na(css$h_ss),]
consensus <- consensus[consensus$h_tx %in% sp_tx,]
consensus$aa <- floor(consensus$h_ss/3) + 1
print(nrow(consensus))
## plot from start
ggplot(consensus, aes(x=aa)) + geom_histogram(binwidth = 1) + xlim(0,100)
ggsave(file = c(file.path(save_plot_path, paste0(org, ".png"))))


