### compare transcripts with stall sites
### how many tx with stall sites have homologs (tx/gene name TXT files)
### read biomart_export and consensus_df

### call stall sites again? without picking the longest before (the conservation script will take care of it)


consensus_yeast <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/consensus_df/yeast.csv"))
consensus_fruitfly <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/consensus_df/fruitfly.csv"))
consensus_zebrafish <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/consensus_df/zebrafish.csv"))
consensus_mouse <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/consensus_df/mouse.csv"))
consensus_human <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/consensus_df/human.csv"))

biomart_yeast <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/biomart_export/yeast.txt"))
biomart_fruitfly <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/biomart_export/fruitfly.txt"))
biomart_zebrafish <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/biomart_export/zebrafish.txt"))
biomart_mouse <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/biomart_export/mouse.txt"))
biomart_human <- read.csv(file.path("/Volumes/USELESS/STALLING/conservation/biomart_export/human.txt"))

biomart_mouse$Gene.name <- toupper(biomart_mouse$Gene.name)

bckg <- toupper(biomart_human$Gene.name)[toupper(biomart_human$Gene.name) %in% toupper(biomart_mouse$Gene.name)]
bckg <- c(bckg, toupper(biomart_human$Gene.name)[toupper(biomart_human$Gene.name) %in% toupper(biomart_zebrafish$Gene.name)])
bckg <- c(bckg, toupper(biomart_human$Gene.name)[toupper(biomart_human$Gene.name) %in% toupper(biomart_fruitfly$Gene.name)])
bckg <- c(bckg, toupper(biomart_human$Gene.name)[toupper(biomart_human$Gene.name) %in% toupper(biomart_yeast$Gene.name)])
bckg <- unique(bckg)

fileConn<-file("/Volumes/USELESS/STALLING/PLOTS/conservation/bckg_human.txt")
writeLines(bckg, fileConn)
close(fileConn)

