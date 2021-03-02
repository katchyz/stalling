### MOTIF p-values

# sample randomly 2397 33nt-long sequences, count motif matches
# (do it multiple times, average?)
fasta_conserved <- read.fasta("../../DATA/sequence/human_conserved.fa")
# subset cds
subset_cds <- function(tx, fasta) {
  utr5_len <- txLengths_h[txLengths_h$tx_name == tx,]$utr5_len
  cds_len <- txLengths_h[txLengths_h$tx_name == tx,]$cds_len
  fa <- fasta[[tx]][(utr5_len+1):(utr5_len+cds_len)]
}

fasta_cds <- sapply(names(fasta_conserved), function(x){subset_cds(x, fasta_conserved)})
# write.fasta(sequences = fasta_cds, names = names(fasta_cds), file.out = "../../DATA/sequence/human_conserved_CDS.fa")
tx_no_ss <- table(substr(names(read.fasta("../../DATA/sequence/CSSs.fa")), 1, 15))

substr_matches <- c()
for (k in 1:1000) {
  substrings <- c()
  for (i in 1:length(tx_no_ss)) {
    tx = names(tx_no_ss[i])
    no_ss = as.numeric(tx_no_ss[i])
    fa = toupper(paste0(fasta_cds[[tx]], collapse = ""))
    for (j in 1:no_ss) {
      # sample string
      start = sample(seq(1, (nchar(fa)-33) ,3), 1)
      substring <- substr(fa, start, start + 32)
      substrings <- c(substrings, substring)
    }
  }
  substr_matches <- c(substr_matches, sum(unlist(gregexpr(motif, substrings)) == 1))
}

mean(substr_matches) # 989

m <- matrix(data = c(1419, 989, 2397-1419, 2397-989), nrow = 2, ncol = 2)
chisq <- chisq.test(m)


# get stall sites, count motif matches
stalls <- as.character(read.csv("../../DATA/sequence/CSSs.txt", sep = "\n", header = F)$V1)
sum(unlist(gregexpr(motif, stalls)) == 1) # 1419 matches out of 2397

# p-values for motif
bases = c("A","C","G","T")
motif = "[ACTG]{15}[GCA][GCA][ACTG]{1}[GCA][GCA][ACTG]{13}"
genomeLen = length(unlist(fasta_conserved))

# motif = "AGGATCTAACCG"
# genomeLen = 1E6

# replicate often enough to get a reasonable precision
n = replicate(1000, {
  # sample a new random genome sequence of the specified size
  genomeSeq = sample(bases, size = genomeLen, replace = TRUE)
  # convert the vector of characters into a character string
  genomeSeq = paste0(genomeSeq, collapse="")
  # find all the occurences of motif within genomeSeq
  matches = gregexpr(motif,genomeSeq)[[1]]
  # return the number of these occurences
  if (matches[1]==-1) 0 else length(matches)
})

# these are the estimated probabilities of 0, 1, 2, ... k occurrences:
table(n)/length(n)

