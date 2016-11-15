ggplot(data, aes(x=X, y=Dome)) + geom_bar(stat="identity") + coord_flip()
ggplot(dome, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()


dupa <- transform(dupa, codon = reorder(codon, value))



shield <- data.frame(data$X, data$Shield)
colnames(shield) <- c("codon", "value")
shield <- transform(shield, codon = reorder(codon, value))
ggplot(shield, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()

##
d <- data.frame(data$X, data$X28hpf)
colnames(d) <- c("codon", "value")
d <- transform(d, codon = reorder(codon, value))
ggplot(d, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()

###
# A-site
data <- read.csv('/Volumes/USELESS/META/codons/codons_A_site.csv', header = TRUE)
d <- data.frame(data$X, data$Shield)
colnames(d) <- c("codon", "value")
d <- transform(d, codon = reorder(codon, value))
ggplot(d, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()


### MOUSE
data <- read.csv('/Volumes/USELESS/META/codons/codons_mouse.csv', header = TRUE)
d <- data.frame(data$X, data$mouse_none)
colnames(d) <- c("codon", "value")
d <- transform(d, codon = reorder(codon, value))
ggplot(d, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()


### Giraldez
data <- read.csv('/Volumes/USELESS/META/codons/codons_Giraldez_A_site.csv', header = TRUE)
d <- data.frame(data$X, data$RPF_12h_run3)
colnames(d) <- c("codon", "value")
d <- transform(d, codon = reorder(codon, value))
ggplot(d, aes(x=codon, y=value)) + geom_bar(stat="identity") + coord_flip()