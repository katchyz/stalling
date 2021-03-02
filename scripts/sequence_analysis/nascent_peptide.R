### positive/negative/basic amino acids before stall site vs before random position

### add at specific positions
### add total

library(seqinr)
library(GenomicFeatures)
library(data.table)
library(ggplot2)

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
gtf_path <- "../../DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("../../DATA/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz",
                                      "../../DATA/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.pep.all.fa.gz",
                                      "../../DATA/fasta/danio_rerio/GRCz10/pep/Danio_rerio.GRCz10.pep.all.fa.gz",
                                      "../../DATA/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz",
                                      "../../DATA/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"))

amino_acids <- list(positive = c("r", "h", "k"),
                    negative = c("d", "e"),
                    uncharged = c("s", "t", "n", "q"),
                    special = c("c", "u", "g", "p"),
                    hydrophobic = c("a", "i", "l", "m", "f", "w", "y", "v"))

plot_path <- "../../sequence"


####### CSSs

gtf_path <- "../../DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]

css <- read.csv(file = "../../DATA/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
css_human <- css[!is.na(css$h_ss),]
css_human$aa <- floor(css_human$h_ss/3)+1
css_human <- css_human[css_human$aa > 30,]
css_human$cds_len <- sapply(css_human$h_tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

fasta_pep <- read.fasta("../../DATA/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz")
names(fasta_pep) <- sapply(getAnnot(fasta_pep), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})
fasta_pep <- fasta_pep[names(fasta_pep) %in% unique(as.character(css_human$h_tx))]

nascent_peptide <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_human$h_tx), css_human$aa))
nascent_peptide <- as.data.table(nascent_peptide)

for(col in names(nascent_peptide)) set(nascent_peptide, i=which(nascent_peptide[[col]] %in% amino_acids$positive), j=col, value="P")
for(col in names(nascent_peptide)) set(nascent_peptide, i=which(nascent_peptide[[col]] %in% amino_acids$negative), j=col, value="N")
for(col in names(nascent_peptide)) set(nascent_peptide, i=which(nascent_peptide[[col]] %in% amino_acids$uncharged), j=col, value="U")
for(col in names(nascent_peptide)) set(nascent_peptide, i=which(nascent_peptide[[col]] %in% amino_acids$special), j=col, value="S")
for(col in names(nascent_peptide)) set(nascent_peptide, i=which(nascent_peptide[[col]] %in% amino_acids$hydrophobic), j=col, value="H")

AAs <- data.table(positive = colSums(nascent_peptide == "P"),
                  negative = colSums(nascent_peptide == "N"),
                  uncharged = colSums(nascent_peptide == "U"),
                  special = colSums(nascent_peptide == "S"),
                  hydrophobic = colSums(nascent_peptide == "H"))

AAs <- AAs / nrow(nascent_peptide) # fraction of sites


### CONTROL
utx <- unique(as.character(css_human$h_tx))
random_pool <- list()
for (i in 1:length(utx)) {
  tx <- utx[i]
  cons <- css_human[css_human$h_tx == tx,]
  cds_len_aa <- floor(cons[1,]$cds_len / 3)
  around_ss <- c()
  for (ss in cons$aa) {
    around_ss <- c(around_ss, c((ss-31):(ss+5)))
  }
  pool <- c(31:(cds_len_aa-5))[!c(31:(cds_len_aa-5)) %in% around_ss]
  random_pool[[tx]] <- pool
}

if (sum(sapply(random_pool, function(x){length(x)}) == 0) > 0) {
  css_human <- css_human[!css_human$h_tx == names(random_pool[sapply(random_pool, function(x){length(x)}) == 0]),]
}


#### control 1000
random1000_pos <- data.table(codon = c(-30:0))
random1000_neg <- data.table(codon = c(-30:0))
random1000_spe <- data.table(codon = c(-30:0))

for (j in 1:10000) {
  css_human$random1 <- sapply(as.character(css_human$h_tx), function(x){sample(random_pool[[x]],1)})
  nascent_peptide_random1 <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_human$h_tx), css_human$random1))
  nascent_peptide_random1 <- as.data.table(nascent_peptide_random1)
  
  for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$positive), j=col, value="P")
  for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$negative), j=col, value="N")
  for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$uncharged), j=col, value="U")
  for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$special), j=col, value="S")
  for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$hydrophobic), j=col, value="H")
  
  AAs_random_positive <- data.table(positive = colSums(nascent_peptide_random1 == "P"))
  AAs_random_negative <- data.table(negative = colSums(nascent_peptide_random1 == "N"))
  AAs_random_special <- data.table(special = colSums(nascent_peptide_random1 == "S"))
  
  AAs_random_positive <- AAs_random_positive / nrow(nascent_peptide_random1)
  AAs_random_negative <- AAs_random_negative / nrow(nascent_peptide_random1)
  AAs_random_special <- AAs_random_special / nrow(nascent_peptide_random1)
  
  random1000_pos <- cbind(random1000_pos, AAs_random_positive$positive)
  random1000_neg <- cbind(random1000_neg, AAs_random_negative$negative)
  random1000_spe <- cbind(random1000_spe, AAs_random_special$special)
  
  # AAs_random1 <- data.table(positive = colSums(nascent_peptide_random1 == "P"),
  #                           negative = colSums(nascent_peptide_random1 == "N"),
  #                           uncharged = colSums(nascent_peptide_random1 == "U"),
  #                           special = colSums(nascent_peptide_random1 == "S"),
  #                           hydrophobic = colSums(nascent_peptide_random1 == "H"))
  
  # AAs_random1 <- AAs_random1 / nrow(nascent_peptide_random1) # fraction of sites
  # AAs_random1$iteration <- rep(j, nrow(AAs_random1))
  print(j)
}

dt <- data.table(codon = c(c(-30:0),c(-30:0)),
                 type = c(rep("CSS", 31), rep("control", 31)),
                 pos = c(AAs$positive, apply(random1000_pos[,2:10001], 1, mean)),
                 neg = c(AAs$negative, apply(random1000_neg[,2:10001], 1, mean)),
                 spe = c(AAs$special, apply(random1000_spe[,2:10001], 1, mean)))
                 

ggplot(dt, aes(x=codon, y=pos, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic()
ggsave(c(file.path(plot_path, "FIXED_CSS_positive.png")))

ggplot(dt, aes(x=codon, y=neg, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic()
ggsave(c(file.path(plot_path, "FIXED_CSS_negative.png")))

ggplot(dt, aes(x=codon, y=spe, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic()
ggsave(c(file.path(plot_path, "FIXED_CSS_special.png")))

# save workspace
# save.image("../../DATA/sequence/nascent_peptide.RData")
# load("../../DATA/sequence/nascent_peptide.RData")

### shaded region with 5th and 95th quantile
dt_pos <- dt
dt_pos$q5 <- rep(apply(random1000_pos[,2:10001], 1, function(x) quantile(x, 0.025)), 2)
dt_pos$q95 <- rep(apply(random1000_pos[,2:10001], 1, function(x) quantile(x, 0.975)), 2)

ggplot(dt_pos, aes(x=codon, y=pos)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_pos, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_CSS_positive.png")))

dt_neg <- dt
dt_neg$q5 <- rep(apply(random1000_neg[,2:10001], 1, function(x) quantile(x, 0.025)), 2)
dt_neg$q95 <- rep(apply(random1000_neg[,2:10001], 1, function(x) quantile(x, 0.975)), 2)

ggplot(dt_neg, aes(x=codon, y=neg)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_neg, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_CSS_negative.png")))

dt_spe <- dt
dt_spe$q5 <- rep(apply(random1000_spe[,2:10001], 1, function(x) quantile(x, 0.025)), 2)
dt_spe$q95 <- rep(apply(random1000_spe[,2:10001], 1, function(x) quantile(x, 0.975)), 2)

ggplot(dt_spe, aes(x=codon, y=spe)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_spe, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_CSS_special.png")))

######

css_human$random1 <- sapply(as.character(css_human$h_tx), function(x){sample(random_pool[[x]],1)})
css_human$random2 <- sapply(as.character(css_human$h_tx), function(x){sample(random_pool[[x]],1)})

# control1
nascent_peptide_random1 <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_human$h_tx), css_human$random1))
nascent_peptide_random1 <- as.data.table(nascent_peptide_random1)

for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$positive), j=col, value="P")
for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$negative), j=col, value="N")
for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$uncharged), j=col, value="U")
for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$special), j=col, value="S")
for(col in names(nascent_peptide_random1)) set(nascent_peptide_random1, i=which(nascent_peptide_random1[[col]] %in% amino_acids$hydrophobic), j=col, value="H")

AAs_random1 <- data.table(positive = colSums(nascent_peptide_random1 == "P"),
                          negative = colSums(nascent_peptide_random1 == "N"),
                          uncharged = colSums(nascent_peptide_random1 == "U"),
                          special = colSums(nascent_peptide_random1 == "S"),
                          hydrophobic = colSums(nascent_peptide_random1 == "H"))

AAs_random1 <- AAs_random1 / nrow(nascent_peptide_random1) # fraction of sites

# control2
nascent_peptide_random2 <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_human$h_tx), css_human$random2))
nascent_peptide_random2 <- as.data.table(nascent_peptide_random2)

for(col in names(nascent_peptide_random2)) set(nascent_peptide_random2, i=which(nascent_peptide_random2[[col]] %in% amino_acids$positive), j=col, value="P")
for(col in names(nascent_peptide_random2)) set(nascent_peptide_random2, i=which(nascent_peptide_random2[[col]] %in% amino_acids$negative), j=col, value="N")
for(col in names(nascent_peptide_random2)) set(nascent_peptide_random2, i=which(nascent_peptide_random2[[col]] %in% amino_acids$uncharged), j=col, value="U")
for(col in names(nascent_peptide_random2)) set(nascent_peptide_random2, i=which(nascent_peptide_random2[[col]] %in% amino_acids$special), j=col, value="S")
for(col in names(nascent_peptide_random2)) set(nascent_peptide_random2, i=which(nascent_peptide_random2[[col]] %in% amino_acids$hydrophobic), j=col, value="H")

AAs_random2 <- data.table(positive = colSums(nascent_peptide_random2 == "P"),
                          negative = colSums(nascent_peptide_random2 == "N"),
                          uncharged = colSums(nascent_peptide_random2 == "U"),
                          special = colSums(nascent_peptide_random2 == "S"),
                          hydrophobic = colSums(nascent_peptide_random2 == "H"))

AAs_random2 <- AAs_random2 / nrow(nascent_peptide_random2) # fraction of sites

AAs$position <- rep("stall_site", nrow(AAs))
AAs$codon <- c(-30:0)
AAs_random1$position <- rep("random1", nrow(AAs_random1))
AAs_random1$codon <- c(-30:0)
AAs_random2$position <- rep("random2", nrow(AAs_random2))
AAs_random2$codon <- c(-30:0)

AAs_all <- rbind(AAs, AAs_random1, AAs_random2)

ggplot(AAs_all, aes(x=codon, y=positive, colour=position)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388"), name = "position on mRNA", labels = c("random 1", "random2", "stall site")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()
ggsave(c(file.path(plot_path, "CSS_positive.png")))

ggplot(AAs_all, aes(x=codon, y=negative, colour=position)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388"), name = "position on mRNA", labels = c("random 1", "random2", "stall site")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()
ggsave(c(file.path(plot_path, "CSS_negative.png")))

ggplot(AAs_all, aes(x=codon, y=uncharged, colour=position)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388"), name = "position on mRNA", labels = c("random 1", "random2", "stall site")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()
ggsave(c(file.path(plot_path, "CSS_uncharged.png")))

ggplot(AAs_all, aes(x=codon, y=special, colour=position)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388"), name = "position on mRNA", labels = c("random 1", "random2", "stall site")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()
ggsave(c(file.path(plot_path, "CSS_special.png")))

ggplot(AAs_all, aes(x=codon, y=hydrophobic, colour=position)) + geom_line(aes(colour = position, linetype = position)) +
  scale_color_manual(values=c("gray", "gray", "#615388"), name = "position on mRNA", labels = c("random 1", "random2", "stall site")) +
  scale_linetype_manual(values = c("dotted", "dotted", "solid")) + theme_classic()
ggsave(c(file.path(plot_path, "CSS_hydrophobic.png")))


################## Ps in all 31 aas
nascent_peptide$Ps <- apply(nascent_peptide, 1, function(x){sum(x == "P")})
nascent_peptide_random1$Ps <- apply(nascent_peptide_random1, 1, function(x){sum(x == "P")})
nascent_peptide_random2$Ps <- apply(nascent_peptide_random2, 1, function(x){sum(x == "P")})

p31 <- data.table(position = c(rep("stall", nrow(nascent_peptide)),
                      rep("random1", nrow(nascent_peptide_random1)),
                      rep("random2", nrow(nascent_peptide_random2))),
           Ps = c(nascent_peptide$Ps, nascent_peptide_random1$Ps, nascent_peptide_random2$Ps))

ggplot(p31, aes(x = Ps, fill = position)) + geom_histogram(position = "dodge")
ggplot(p31, aes(x = Ps, fill = position)) + geom_density(alpha = 0.3)

s <- p31[p31$position == "stall",]
ggplot(s, aes(x = Ps)) + geom_histogram(binwidth = 1)
r1 <- p31[p31$position == "random1",]
ggplot(r1, aes(x = Ps)) + geom_histogram(binwidth = 1)
r2 <- p31[p31$position == "random2",]
ggplot(r2, aes(x = Ps)) + geom_histogram(binwidth = 1)

### Ps in the last 10 AAs
nascent_peptide$Ps11 <- apply(nascent_peptide[,21:31], 1, function(x){sum(x == "P")})
nascent_peptide_random1$Ps11 <- apply(nascent_peptide_random1[,21:31], 1, function(x){sum(x == "P")})
nascent_peptide_random2$Ps11 <- apply(nascent_peptide_random2[,21:31], 1, function(x){sum(x == "P")})

p11 <- data.table(position = c(rep("stall", nrow(nascent_peptide)),
                               rep("random1", nrow(nascent_peptide_random1)),
                               rep("random2", nrow(nascent_peptide_random2))),
                  Ps = c(nascent_peptide$Ps11, nascent_peptide_random1$Ps11, nascent_peptide_random2$Ps11))

ggplot(p11, aes(x = Ps, fill = position)) + geom_histogram(position = "dodge")
s <- p11[p11$position == "stall",]
ggplot(s, aes(x = Ps)) + geom_histogram(binwidth = 1)
r1 <- p11[p11$position == "random1",]
ggplot(r1, aes(x = Ps)) + geom_histogram(binwidth = 1)
r2 <- p11[p11$position == "random2",]
ggplot(r2, aes(x = Ps)) + geom_histogram(binwidth = 1)



######################
## export those without p/g/d at P-site, e at A-site, d/e at -2

nascent_peptide <- t(mapply(function(x,y){fasta_pep[[x]][(y-15):(y+15)]}, as.character(css_human$h_tx), css_human$aa))
nascent_peptide <- as.data.table(nascent_peptide)
rownames(nascent_peptide) <- paste0(css_human$h_tx, "_", css_human$aa)

REST <- rownames(nascent_peptide)[nascent_peptide$V16 != "p" & nascent_peptide$V16 != "g" & nascent_peptide$V16 != "d" &
                  nascent_peptide$V17 != "e" & nascent_peptide$V14 != "d" & nascent_peptide$V14 != "e"]
# the rest!
save(REST, file = "../../DATA/structure/rest_not_explained_by_AA.Rsave")

nascent_peptide[(nascent_peptide$V14 == "e" | nascent_peptide$V14 == "d") &
                  (nascent_peptide$V16 == "p" | nascent_peptide$V16 == "g" | nascent_peptide$V16 == "d")]
# 227 cases -2 and P-site

nascent_peptide[(nascent_peptide$V17 == "e") &
                  (nascent_peptide$V16 == "p" | nascent_peptide$V16 == "g" | nascent_peptide$V16 == "d")]
# 125 e at A-site and either p/g/d at P-site
nascent_peptide[(nascent_peptide$V17 == "e") &
                  (nascent_peptide$V16 != "p" | nascent_peptide$V16 != "g" | nascent_peptide$V16 != "d")]
# 326 alone e at A-site



