## export consensus to python
library(GenomicFeatures)

#consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/median/8"
#consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/threshold/8"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
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


for (i in 1:length(organisms)) {
  org <- organisms[i]
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus$ss <- consensus$ss - 1 # 0-based coord
  
  #write.csv(consensus, file.path("/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median", paste0(org, ".csv")))
  write.csv(consensus, file.path("/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median_all_peaks", paste0(org, ".csv")))
}



