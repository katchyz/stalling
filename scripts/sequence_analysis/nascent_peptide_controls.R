### nascent peptide (pos/neg/special) analysis for subsets (cons3, y, f, z, m)


load(file = "../../DATA/sequence/cons3_human.Rsave") # c3h
css_h <- css[!is.na(css$h_ss),]
css_h$aa <- floor(css_h$h_ss/3)+1
css_human <- css_human[css_human$aa > 30,]

# human, in 3 organisms
c3h <- c3h[c3h$ss_aa > 30]
c3h$cds_len <- sapply(c3h$tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

fasta_pep <- read.fasta("../../DATA/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz")
names(fasta_pep) <- sapply(getAnnot(fasta_pep), function(x){strsplit(substring(x, regexpr("transcript:", x) + 11), split = " ")[[1]][1]})
fasta_pep <- fasta_pep[names(fasta_pep) %in% unique(as.character(c3h$tx))]

# mouse
css_m <- css[!is.na(css$m_ss),]
css_m$aa <- floor(css_m$m_ss/3)+1
css_m <- css_m[css_m$aa > 30]
css_m$cds_len <- sapply(css_m$m_tx, function(x){txLengths_m[txLengths_m$tx_name == x,]$cds_len})

fasta_pep <- fasta_mouse[names(fasta_mouse) %in% unique(as.character(css_m$m_tx))]

# zebrafish
css_z <- css[!is.na(css$z_ss),]
css_z$aa <- floor(css_z$z_ss/3)+1
css_z <- css_z[css_z$aa > 30]
css_z$cds_len <- sapply(css_z$z_tx, function(x){txLengths_z[txLengths_z$tx_name == x,]$cds_len})

fasta_pep <- fasta_zebra[names(fasta_zebra) %in% unique(as.character(css_z$z_tx))]

# fruitfly
css_f <- css[!is.na(css$f_ss),]
css_f$aa <- floor(css_f$f_ss/3)+1
css_f <- css_f[css_f$aa > 30]
css_f$cds_len <- sapply(css_f$f_tx, function(x){txLengths_f[txLengths_f$tx_name == x,]$cds_len})

fasta_pep <- fasta_fruit[names(fasta_fruit) %in% unique(as.character(css_f$f_tx))]

# yeast
css_y <- css[!is.na(css$y_ss),]
css_y$aa <- floor(css_y$y_ss/3)+1
css_y <- css_y[css_y$aa > 30]
css_y$cds_len <- sapply(css_y$y_tx, function(x){txLengths_y[txLengths_y$tx_name == x,]$cds_len})

fasta_pep <- fasta_yeast[names(fasta_yeast) %in% unique(as.character(css_y$y_tx))]

###
nascent_peptide <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_y$y_tx), css_y$aa))
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
utx <- unique(as.character(css_y$y_tx))
random_pool <- list()
for (i in 1:length(utx)) {
  t <- utx[i]
  cons <- css_y[css_y$y_tx == t,]
  cds_len_aa <- floor(cons[1,]$cds_len / 3)
  around_ss <- c()
  for (ss in cons$aa) {
    around_ss <- c(around_ss, c((ss-31):(ss+5)))
  }
  pool <- c(31:(cds_len_aa-5))[!c(31:(cds_len_aa-5)) %in% around_ss]
  random_pool[[t]] <- pool
}

if (sum(sapply(random_pool, function(x){length(x)}) == 0) > 0) {
  css_y <- css_y[!css_y$y_tx == names(random_pool[sapply(random_pool, function(x){length(x)}) == 0]),]
}


#### control 1000
random1000_pos <- data.table(codon = c(-30:0))
random1000_neg <- data.table(codon = c(-30:0))
random1000_spe <- data.table(codon = c(-30:0))

for (j in 1:1000) {
  css_y$random1 <- sapply(as.character(css_y$y_tx), function(x){sample(random_pool[[x]],1)})
  nascent_peptide_random1 <- t(mapply(function(x,y){fasta_pep[[x]][(y-30):y]}, as.character(css_y$y_tx), css_y$random1))
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
                 pos = c(AAs$positive, apply(random1000_pos[,2:1001], 1, mean)),
                 neg = c(AAs$negative, apply(random1000_neg[,2:1001], 1, mean)),
                 spe = c(AAs$special, apply(random1000_spe[,2:1001], 1, mean)))




plot_path <- "../../DATA/sequence"
### shaded region with 5th and 95th quantile
dt_pos <- dt
dt_pos$q5 <- rep(apply(random1000_pos[,2:1001], 1, function(x) quantile(x, 0.025, na.rm = T)), 2)
dt_pos$q95 <- rep(apply(random1000_pos[,2:1001], 1, function(x) quantile(x, 0.975, na.rm = T)), 2)

ggplot(dt_pos, aes(x=codon, y=pos)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_pos, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_YEAST_positive.png")))

dt_neg <- dt
dt_neg$q5 <- rep(apply(random1000_neg[,2:1001], 1, function(x) quantile(x, 0.025)), 2)
dt_neg$q95 <- rep(apply(random1000_neg[,2:1001], 1, function(x) quantile(x, 0.975)), 2)

ggplot(dt_neg, aes(x=codon, y=neg)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_neg, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_YEAST_negative.png")))

dt_spe <- dt
dt_spe$q5 <- rep(apply(random1000_spe[,2:1001], 1, function(x) quantile(x, 0.025)), 2)
dt_spe$q95 <- rep(apply(random1000_spe[,2:1001], 1, function(x) quantile(x, 0.975)), 2)

ggplot(dt_spe, aes(x=codon, y=spe)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), name = "type") +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_ribbon(data=subset(dt_spe, type == "CSS"), aes(ymin=q5, ymax=q95),fill="black", alpha=0.2)
ggsave(c(file.path(plot_path, "SHADE_YEAST_special.png")))
