### CONSERVATION
library(data.table)

## get list of homologs
## get list of conserved stall sites

## (align transcripts of interest)
## get position on CDS (codon number)

homologs <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/homologs.csv", header = F)
colnames(homologs) <- c("gene", "yeast", "fruitfly", "zebrafish", "mouse", "human")
setDT(homologs)

conserved <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
colnames(conserved) <- c("gene", "yeast_tx", "yeast_pos", "fruitfly_tx", "fruitfly_pos", "zebrafish_tx", "zebrafish_pos",
                         "mouse_tx", "mouse_pos", "human_tx", "human_pos")
setDT(conserved)
conserved$gene <- sapply(strsplit(as.character(conserved$gene), split = "_"), function(x){x[1]})

### get all homologs (with or without stall site) in one table

conserved$yeast_tx <- sapply(conserved$gene, function(x){homologs[homologs$gene == x,]$yeast})
conserved$fruitfly_tx <- sapply(conserved$gene, function(x){homologs[homologs$gene == x,]$fruitfly})
conserved$zebrafish_tx <- sapply(conserved$gene, function(x){homologs[homologs$gene == x,]$zebrafish})
conserved$mouse_tx <- sapply(conserved$gene, function(x){homologs[homologs$gene == x,]$mouse})
conserved$human_tx <- sapply(conserved$gene, function(x){homologs[homologs$gene == x,]$human})

### pos are in CDS-coordinates, 0-based
### to convert to codons: pos / 3 + 1

conserved$yeast_pos <- floor(conserved$yeast_pos / 3 + 1)
conserved$fruitfly_pos <- floor(conserved$fruitfly_pos / 3 + 1)
conserved$zebrafish_pos <- floor(conserved$zebrafish_pos / 3 + 1)
conserved$mouse_pos <- floor(conserved$mouse_pos / 3 + 1)
conserved$human_pos <- floor(conserved$human_pos / 3 + 1)

colnames(conserved) <- c("gene", "yeast_tx", "yeast_codon", "fruitfly_tx", "fruitfly_codon", "zebrafish_tx", "zebrafish_codon",
                         "mouse_tx", "mouse_codon", "human_tx", "human_codon")

write.csv(conserved, "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_homologs.csv")


