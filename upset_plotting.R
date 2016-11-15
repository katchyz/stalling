upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"), list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"), list(plot=histogram, x="ReleaseDate")), ncols = 2))


upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")

library(UpSetR)
data <- read.csv('/Volumes/USELESS/META/CONS_gene_peak.csv', header = TRUE)
colnames(data) <- c("X", "zebrafish", "yeast", "mouse", "zebrafish II", "fruit fly")

upset(data, sets = c("our", "giraldez", "mouse", "fruitfly", "yeast"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")

data$giraldez <- NULL

upset(data, sets = c("zebrafish", "mouse", "fruit fly", "yeast"), sets.bar.color = "#D55853", order.by = "freq", name.size = 16, matrix.color = "#D55853")