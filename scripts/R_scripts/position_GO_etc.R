# load consensus
# subtract utr5_len
# plot distribution of peaks
library(ggplot2)

consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
#consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median_all_peaks"
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/start"
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

## start
for (i in 1:length(organisms)) {
  org <- organisms[i]
  ## get utr5 len
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  txLengths <- txLengths[order(txLengths$tx_name),]
  ## get consensus, add utr5_len
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) ## consensus
  consensus$utr5_len <- sapply(consensus$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
  consensus$ss <- consensus$start - consensus$utr5_len
  # exclude first 15
  consensus <- consensus[consensus$ss > 15,]
  save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))
  print(org)
  print(nrow(consensus))
  ## plot, save
  ggplot(consensus, aes(x=ss)) + geom_histogram() #+ scale_x_log10()
  #ggsave(file = c(file.path(save_plot_path, paste0(org, ".png"))))
}

## distance from stop codon
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/signal_peptide/stop"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  ## get utr5 len
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  txLengths <- txLengths[order(txLengths$tx_name),]
  ## get consensus, add utr5_len
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) ## consensus
  consensus$from_stop <- consensus$cds_len - consensus$ss
  print(org)
  ## plot, save
  ggplot(consensus, aes(x=log10(from_stop))) + geom_histogram() + scale_x_reverse()
  ggsave(file = c(file.path(save_plot_path, paste0(org, ".png"))))
}



##########
## get CDS to unit length
# get: cds_len and ss
# plot stall sites at positions ss/cds_len (on a scale 0-1)

for (i in 1:length(organisms)) {
  org <- organisms[i]
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) ## consensus
  consensus$ss_unit <- consensus$ss / consensus$cds_len
  # plot
  ggplot(consensus, aes(x = ss_unit)) + geom_histogram(bins = 50)
  ggsave(file = c(file.path(save_plot_path, "unit", paste0(org, ".png"))))
}

######
go_path <- "/Volumes/USELESS/STALLING/GO"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  txLengths <- txLengths[order(txLengths$tx_name),]
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  txlen <- txlen[order(rownames(txlen)),]
  ##
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) ## consensus
  gene_set <- unique(txlen[txlen$tx_name %in% consensus$seqnames,]$gene_id)
  ##
  fileConn<-file(c(file.path(go_path, paste0(org, "_ss.txt"))))
  writeLines(gene_set, fileConn)
  close(fileConn)
  ##
  bckg <- unique(txlen[!(txlen$tx_name %in% consensus$seqnames),]$gene_id)
  ##
  fileConn<-file(c(file.path(go_path, paste0(org, "_bckg.txt"))))
  writeLines(gene_set, fileConn)
  close(fileConn)
}











############################
############################
###### LOGOS ###############

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))


save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/3lib"
### logos
for (i in 1:length(organisms)) {
  org <- organisms[i]
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus <- consensus[consensus$freq > 2,]
  # get fasta file around these peaks
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-15):(y+17)]}, as.character(consensus$seqnames), consensus$start))
  logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
  logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
  # write to file
  fileConn<-file(c(file.path(save_logo_path, paste0(org, ".txt"))))
  writeLines(logo, fileConn)
  close(fileConn)
}

save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/long"
### logos
for (i in 1:length(organisms)) {
  org <- organisms[i]
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus <- consensus[consensus$start > 90,]
  # get fasta file around these peaks
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-90):(y+17)]}, as.character(consensus$seqnames), consensus$start))
  logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
  logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
  # write to file
  fileConn<-file(c(file.path(save_logo_path, paste0(org, ".txt"))))
  writeLines(logo, fileConn)
  close(fileConn)
}

#### ingolia 2011 - get common in mouse_chx and mouse_none, save logo

