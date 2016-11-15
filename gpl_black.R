data <- read.csv('../meta_data/split_by_length/stall_csv/5dpf.csv', header = TRUE)
data$X <- NULL
data$scale <- c(-15:17)

d <- melt(data, id.var="scale")
ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + scale_fill_brewer(name=element_blank(), palette="RdGy") + theme_black() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
# position="stack"

ggsave(file = "../plots/stall/stall_OTPG/fill/black/7_5dpf.png")

ggplot(f, aes(x=stages, y=value, fill=variable)) + geom_bar(stat="identity", position="dodge") + scale_fill_brewer(name=element_blank(), palette="RdGy") + theme_black() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())