## MSE
# human
library(GenomicRanges)
library(seqinr)

# consensus from 2 libraries
human_2libs <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median/human.csv")
human_2libs <- human_2libs[human_2libs$ss >= 15,]

# conserved in human and others
human_conse <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
colnames(human_conse) <- c("gene", "Y", "Y_ss", "F", "F_ss", "Z", "Z_ss", "M", "M_ss", "H", "H_ss")

human_conse <- human_conse[(!is.na(human_conse$H_ss)) & ( (!is.na(human_conse$Y_ss)) | (!is.na(human_conse$F_ss)) |
                                                            (!is.na(human_conse$Z_ss)) | (!is.na(human_conse$M_ss)) ),]

save(human_2libs, file = "/Users/kasia/Desktop/human_2libs.Rsave")
save(human_conse, file = "/Users/kasia/Desktop/human_conserved.Rsave")

##########
fasta <- "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
fasta_cdna <- read.fasta(fasta)
human_2libs <- human_2libs[human_2libs$start >= 18,]
logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-17):(y+15)]}, as.character(human_2libs$seqnames), human_2libs$start))
save(logo, file = "/Users/kasia/Desktop/human_logo.Rsave")

# mapFromTranscripts

