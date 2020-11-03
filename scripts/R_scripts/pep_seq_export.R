### export pep seq


library(seqinr)
library(GenomicFeatures)
library(data.table)
library(ggplot2)


fasta_pep <- read.fasta("/Volumes/USELESS/DATA/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz")
names(fasta_pep) <- sapply(getAnnot(fasta_pep), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})

pep100 <- human_ss[human_ss$aa > 100,]
fasta_pep <- fasta_pep[names(fasta_pep) %in% pep100$h_tx]

pep_seq <- t(mapply(function(x,y){fasta_pep[[x]][(y-100):(y-1)]}, as.character(pep100$h_tx), pep100$aa))
pep_seq <- as.data.table(pep_seq)
pep_seqs <- do.call(paste0, pep_seq[,1:100])
pep_seqs <- toupper(pep_seqs[nchar(pep_seqs) == 100])

pep100$id <- paste0(pep100$h_tx, "_", pep100$h_ss)

write.fasta(as.list(pep_seqs), pep100$id,
            file.out = "/Volumes/USELESS/STALLING/PLOTS/logo_stall/CSSs_pep.fa")

### export cds seq
nt_seq <- t(mapply(function(x,y){fasta_cds[[x]][(y-15):(y+17)]}, as.character(human_ss$tx), human_ss$nt_cds))
nt_seq <- as.data.table(nt_seq)
css_seqs <- do.call(paste0, nt_seq[,1:33])
css_seqs <- toupper(css_seqs[nchar(css_seqs) == 33])

fileConn<-file("/Volumes/USELESS/STALLING/PLOTS/logo_stall/CSSs.txt")
writeLines(css_seqs, fileConn)
close(fileConn)

write.fasta(as.list(css_seqs), human_ss$tx_ss,
            file.out = "/Volumes/USELESS/STALLING/PLOTS/logo_stall/CSSs.fa")
