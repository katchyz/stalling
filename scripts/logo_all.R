### logos for every lib
library(ggplot2)
library(reshape2)
library(scales)

library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(RColorBrewer)
library(colorRamps)

library(seqinr)

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))

peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/all"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  # load peaks
  libs <- list.files(path = peaks_path, pattern = paste0("^", org))
  
  for (j in 1:length(libs)) {
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- peaks[start > 15]
    
    ### logos - 5
    logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-15):(y+17)]}, as.character(peaks$seqnames), peaks$start))
    logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
    logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
    # write to file
    fileConn<-file(c(file.path(save_logo_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".txt"))))
    writeLines(logo, fileConn)
    close(fileConn)
  }
}

