library(ggplot2)
library(plyr)
library(dplyr)
matrix <- read.csv("/Users/kasia/Documents/PhD/matrix_go_fq.csv", header = TRUE)

css = c("cct3", "rps6", "sf3b4", "mcmbp", "ndnl2", "mcm4", "gtpbp4", "vcp", "snrpb", "snrpa", "arpc2", "aurka", "ddit4", "rpl18a", "gmps", "ilf2", "eif5", "vdac2", "csnk2b", "eif3i", "serinc1", "prkrip1", "suds3")

matrix = within(matrix, {
  cons_ss = ifelse(gene_name %in% css, "C", "N")
})

matrix <- arrange(matrix, cons_ss)
matrix <- arrange(matrix, -row_number())

ggplot(matrix, aes(x=log2(X04_dome_gene_rna), y=log2(X04_dome_gene_ribo))) + geom_point(shape=1) + geom_point(aes(colour=cons_ss),size=3) + geom_text(aes(label=ifelse(cons_ss=="C",as.character(gene_name),'')),hjust=0,just=0, size=4)

ggsave(file = "/Volumes/USELESS/META/CSS.png")