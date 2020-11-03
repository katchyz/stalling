upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"), list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"), list(plot=histogram, x="ReleaseDate")), ncols = 2))


upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")

library(UpSetR)
data <- read.csv('/Volumes/USELESS/META/all_the_crap/CONS_gene_peak.csv', header = TRUE)
colnames(data) <- c("X", "zebrafish", "yeast", "mouse", "zebrafish II", "fruit fly")

upset(data, sets = c("our", "giraldez", "mouse", "fruitfly", "yeast"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")

data$giraldez <- NULL

upset(data, sets = c("zebrafish", "mouse", "fruit fly", "yeast"), sets.bar.color = "#D55853", order.by = "freq", name.size = 16, matrix.color = "#D55853")

### stalling ###
library(UpSetR)

data <- read.csv('/Volumes/USELESS/STALLING/conservation/upsetr.csv', header = TRUE)
data <- read.csv('/Volumes/USELESS/STALLING/conservation/upsetr_median5.csv', header = TRUE)
data <- read.csv('/Volumes/USELESS/STALLING/conservation/upsetr_threshold5.csv', header = TRUE)
data <- read.csv('/Volumes/USELESS/STALLING/conservation/upsetr_median8.csv', header = TRUE)
data <- read.csv('/Volumes/USELESS/STALLING/conservation/upsetr_threshold8.csv', header = TRUE)
data2 <- read.csv('/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median/NEW_median.csv', header = TRUE)
data1 <- read.csv('/Volumes/USELESS/STALLING/conservation/consensus_df_all/NEW_median_all_peaks/NEW_median_all_peaks.csv', header = TRUE)

# upset(data)


upset(data1, order.by = "freq", text.scale = 2, point.size = 5, sets.bar.color = "#615388",
      mainbar.y.label = "No. conserved stall sites", sets.x.label = "No. stall sites")
upset(data2, order.by = "freq", text.scale = 2, point.size = 5, sets.bar.color = "#615388",
      mainbar.y.label = "No. conserved stall sites", sets.x.label = "No. stall sites")

genes <- as.character(data1$gene_ss)
genes <- unique(sapply(strsplit(genes, "_"), function(x){x[1]}))
fileConn<-file("/Volumes/USELESS/STALLING/PLOTS/conservation/genes_with_cons_ss_all.txt")
writeLines(genes, fileConn)
close(fileConn)

genes2 <- as.character(data2$gene_ss)
genes2 <- unique(sapply(strsplit(genes2, "_"), function(x){x[1]}))
fileConn<-file("/Volumes/USELESS/STALLING/PLOTS/conservation/genes_with_cons_ss.txt")
writeLines(genes2, fileConn)
close(fileConn)
