# compare AA frequencies in the whole 40aa upstream of CSSs vs control
library(Biostrings)
library(data.table)
library(ggplot2)
library(seqinr)

#load("~/Downloads/kmer.Rsave")
load("../../DATA/sequence/fasta_pep.Rsave")
load("../../DATA/sequence/CSSs.Rsave")
setDT(CSSs)
CSSs <- CSSs[order(CSSs$tx),]

# get first 33 aa in css and control
#AA_css <- translate(DNAStringSet(sapply(kmer$css_seq, function(x){substr(x,3,101)})))

pep_len <- data.table(tx = names(fasta_pep), pep_len = sapply(fasta_pep, length))
CSSs$pep_len <- sapply(as.character(CSSs$tx), function(x){pep_len[tx == x]$pep_len})

# control positions
ctrl_coord <- list()
for (txn in names(fasta_pep)) {
  #print(txn)
  tx_len <- CSSs[tx == txn]$pep_len[1]
  #print(tx_len)
  ss <- CSSs[tx == txn]$ss_aa
  #print(ss)
  coord_to_excl <- as.numeric(sapply(ss, function(x){c((x-40):(x+40))}))
  coord_to_excl <- coord_to_excl[coord_to_excl > 0]
  #print(coord_to_excl)
  ctrl <- c(1:(tx_len-1))[-coord_to_excl]
  ctrl_coord[[txn]] <- ctrl
}

ctrl_coord <- lapply(ctrl_coord, function(x){x[x > 40]})

tx_to_excl <- names(ctrl_coord[sapply(ctrl_coord, length) == 0])
CSSs <- CSSs[!(tx %in% tx_to_excl)]

CSSs$ctrl <- sapply(as.character(CSSs$tx), function(x){sample(ctrl_coord[[x]], 1)})

# exclude SS closer than 40aa to start
CSSs <- CSSs[ss_aa > 40] # 1802 ss

seq_css <- as.character(mapply(function(x,y){fasta_pep[[x]][(y-40):(y-1)]}, as.character(CSSs$tx), CSSs$ss_aa))
seq_control <- as.character(mapply(function(x,y){fasta_pep[[x]][(y-40):(y-1)]}, as.character(CSSs$tx), CSSs$ctrl))
seq_all <- as.character(unlist(fasta_pep[names(fasta_pep) %in% as.character(CSSs$tx)]))

t_css <- table(seq_css) / length(seq_css)
t_ctr <- table(seq_control) / length(seq_control)
t_all <- table(seq_all) / length(seq_all)
t_all <- t_all[-18] # remove 'u'
t_all <- t_all / sum(t_all)

# genome
fasta <- read.fasta("../../DATA/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz")
fasta_annot <- sapply(getAnnot(fasta), function(x){sub(".*transcript:", "", x)})
fasta_annot <- as.character(sapply(fasta_annot, function(x){sub(" gene_biotype.*", "", x)}))
names(fasta) <- fasta_annot

t_gen <- table(as.character(unlist(fasta))) / length(unlist(fasta))
t_gen <- t_gen[2:23]
t_gen <- t_gen[-21] #remove x
t_gen <- t_gen[-18] #remove u
t_gen <- t_gen / sum(t_gen)

### compare frequencies
res_css <- chisq.test(t_css*1000, p = t_all)
res_control <- chisq.test(t_ctr*1000, p = t_all)

dt = data.table(aa = rep(toupper(names(t_css)), 3), freq = c(t_css, t_ctr, t_all),
                type = c(rep("CSS", 20), rep("control", 20), rep("all", 20)))

ggplot(dt, aes(x = aa, y = freq, fill = type)) + geom_bar(stat = 'identity', position = 'dodge') + theme_classic() +
  scale_fill_manual(values=c("#615388", "gray", "black"), name = " ", labels = c("CSS", "control", "genome"))

# ggsave("nascent_pep_AA_distribution.png", height = 6, width = 12)

# aromatic
aromatic <- c("f", "w", "y")
ar_css <- c(sum(t_css[names(t_css) %in% aromatic]), sum(t_css[!names(t_css) %in% aromatic]))
ar_ctr <- c(sum(t_ctr[names(t_ctr) %in% aromatic]), sum(t_ctr[!names(t_ctr) %in% aromatic]))
ar_all <- c(sum(t_all[names(t_all) %in% aromatic]), sum(t_all[!names(t_all) %in% aromatic]))
ar_gen <- c(sum(t_gen[names(t_gen) %in% aromatic]), sum(t_gen[!names(t_gen) %in% aromatic]))

chisq.test(ar_css*1000, p = ar_gen)
chisq.test(ar_ctr*1000, p = ar_gen)
