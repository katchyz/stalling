setwd("~/Documents/PhD/scripts")
library(reshape2)
library(ggplot2)
# library(RColorBrewer)

data <- read.csv('/Volumes/USELESS/meta_data_MAC/split_by_length/stop_csv/Dome.csv', header = TRUE)
data$X <- NULL
data$scale <- c(-15:17)
colnames(data) <- c("25", "26", "27", "28", "29", "30", "31", "scale")

d <- melt(data, id.var="scale")
ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + scale_fill_brewer(name=element_blank(), palette="RdGy") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(fill=guide_legend(title="footprint\nlength"))

# position="stack"

ggsave(file = "../plots/stop/stop_OTPG/stack/3_Dome.png")

ggsave(file = "/Users/kasia/Documents/PhD/misc./EMBO/fld.png")
