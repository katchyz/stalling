### match GAxGAxGAx pattern
# select well-expressed transcripts
# Ribo-Seq
# Shape-Seq
# GO terms
# is this conserved?
library(seqinr)
library(stringr)
library(GenomicFeatures)
library(ggplot2)

### Zebrafish, GRCz10
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
#fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- lapply(fasta_cds, function(x){paste(x, collapse = "")})

# GA.......GA.GA.GA.GA.G
# scoring against HMM profile ??
GAmatch <- lapply(fasta_cds, function(x){str_locate_all(x, 'ga.......ga.ga.ga.ga.g')})

#sum(sapply(GAmatch, function(x){length(x[[1]]) > 0})) # 3063 tx with matching pattern
GAmatch <- GAmatch[sapply(GAmatch, function(x){length(x[[1]]) > 0})]

# filter those genes by FPKM (GRCz10)
#load("/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load("/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

high_fpkm <- rownames(FPKM_256[FPKM_256$exons_RNA_fpkm > 1,])
GAmatch <- GAmatch[names(GAmatch) %in% high_fpkm]

# single transcript per gene??

### Ribo-Seq (is that GRCz10 ?????)
#load("/Volumes/USELESS/META/SHAPES/riboseq_2_4.Rdata")
load("/Volumes/USELESS/META/SHAPES/riboseq_256.Rdata")
ribo256 <- sapply(riboseq_256, function(x){unname(tapply(x, (seq_along(x)-1) %/% 3, sum))})
ribo256zs <- sapply(ribo256, function(x){c(0, as.numeric(scale(x[2:(length(x)-1)])), 0)})

### Shape-Seq
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
shape256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
shape256 <- split(shape256, shape256$trnames)
names(shape256) <- sapply(names(shape256), function(x){substr(x,1,18)})
shape256 <- shape256[sapply(shape256, function(x){!is.na(sum(x$dtcr))})]
shape256 <- shape256[sapply(shape256, function(x){sum(x$dtcr) > 0})]
shape256 <- sapply(shape256, function(x){x$dtcr/sum(x$dtcr)})

# index in zs: ceiling(x/3)
zs <- list()
shape <- list()
for (n in names(GAmatch)) {
  zs[[n]] <- list()
  shape[[n]] <- list()
  for (j in seq(1,nrow(GAmatch[[n]][[1]]))) {
    start <- GAmatch[[n]][[1]][j,][1]
    end <- GAmatch[[n]][[1]][j,][2]
    start3 <- ceiling(start/3)
    end3 <- ceiling(end/3)
    #print(max(ribo256zs[[n]][start:end])) # on some of them, they're high
    #zs <- append(zs, max(ribo256zs[[n]][start:end]))
    zs[[n]][[j]] <- ribo256zs[[n]][start3:end3]
    shape[[n]][[j]] <- shape256[[n]][start:end]
  }
}

# 15% (almost 300 tx) have z-scores higher than 5; 25% higher than 3

zs <- zs[sapply(zs, function(x){length(x) > 0})]

ss <- list()
for (n in names(zs)) {
  for (i in seq(1, length(zs[[n]]))) {
    if (sum(is.na(zs[[n]][[i]])) == 0) {
      if (max(zs[[n]][[i]][3:4]) > 3) {
        print(n)
        print(i)
        ss[[n]] <- i
      }
    }
  }
}

######
# plot shape (sum meta data)
shape <- shape[sapply(shape, function(x){length(x) > 0})]

## for SS
sum_shape_ss <- rep(0,22)
for (n in names(ss)) {
  if (length(shape[[n]][[ss[[n]]]]) == 22) {
    sum_shape_ss <- sum_shape_ss + shape[[n]][[ss[[n]]]]
  }
}

sum_shape_ss <- data.frame(sum_shape_ss, scale = seq(1,22))
ggplot(sum_shape_ss, aes(x=scale, y=sum_shape_ss)) + geom_bar(stat = "identity")
ggsave(file = "/Users/kasia/Desktop/shape_GA_ss.png")

## for all with GA-pattern
shape_noss <- shape[!names(shape) %in% names(ss)]

sum_shape <- rep(0,22)
#for (n in sample(names(shape_noss), 100)) {
for (n in names(shape_noss)) {
  for (i in seq(1, length(shape_noss[[n]]))) {
    if (length(shape_noss[[n]][[i]]) == 22) {
      if (sum(is.na(shape_noss[[n]][[i]])) == 0) {
        sum_shape <- sum_shape + shape_noss[[n]][[i]]
      }
    }
  }
}

sum_shape <- data.frame(sum_shape, scale = seq(1,22))
ggplot(sum_shape, aes(x=scale, y=sum_shape)) + geom_bar(stat = "identity")
ggsave(file = "/Users/kasia/Desktop/shape_GA_noss.png")


# split by: SS (is there or is not) - are the profiles any different?


