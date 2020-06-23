### STALLING - analyze human libs
library(seqinr)
library(GenomicFeatures)
library(dplyr)
library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(colorRamps)

# fasta (add to data table)
#fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz")

# get longest tx
txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# get transcript with longest cds per gene
txlen <- arrange(txLengths, gene_id, desc(cds_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]
# cds lengths multiple of 3
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

## for subsetting CDS by txlen coordinates 
txlen <- txlen[txlen$cds_len %%3 == 0,] # get the cds which length is multiple of 3
txlen <- txlen[rownames(txlen) %in% names(fasta_cds),] # get the tx that are in fasta
fasta_cds <- fasta_cds[names(fasta_cds) %in% rownames(txlen)]
fasta_cds <- fasta_cds[order(names(fasta_cds))]
fasta_cds <- fasta_cds[sapply(fasta_cds, function(x){length(x)}) == txlen$cds_len] ## check if fasta and txlen are equal
txlen <- txlen[rownames(txlen) %in% names(fasta_cds),]

#cds_coord <- mapply(function(x,y){(x+1):(x+y)}, txlen$utr5_len, txlen$cds_len)
### cds_coord without first 15 and last 5 codons (Ingolia et al., 2011)
cds_coord <- mapply(function(x,y){(x+46):(x+y-15)}, txlen$utr5_len, txlen$cds_len)
tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
grdf_coord <- unlist(grdf_coord)

#######################################################################################################
libs_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges"
peaks_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/peaks_excl_15_5"
libs <- list.files(path = libs_path)

for (i in 1:length(libs)) {
  ## load file
  load(file = c(file.path(libs_path, libs[i])))
  #load("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges/andreev2015_SRR1173905.Rsave") # gr
  riboseq_lengths <- colnames(mcols(gr))[2:length(colnames(mcols(gr)))]
  grdf <- as.data.frame(gr)
  setDT(grdf)
  
  rm(gr)
  
  grdf <- grdf[seqnames %in% rownames(txlen),] # get the longest CDSs, multiple of 3 only
  grdf <- grdf[grdf_coord,] # get CDS coordinates
  
  ## add fasta... 
  #grdf$nt <- as.character(unlist(fasta_cds))
  
  ## aggregate by codons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #df <- data.frame(seq=c('a','a','a','a','a','a','b','b','b'), riboseq = c(1,2,3,4,5,6,4,3,2), dupa = c(2,3,4,6,7,8,1,2,3), nt = c('a', 't', 'g', 'c', 'a', 'c', 'a', 't', 'g'))
  #setDT(df)
  #df[, .(riboseq = sum(riboseq), dupa = sum(dupa), codon = paste(nt, sep="", collapse="")), by= (seq(nrow(df)) - 1) %/% 3]
  #cols <- paste0("X", riboseq_lengths)
  #df[, lapply(.SD, sum), by = (seq(nrow(df)) - 1) %/% 3, .SDcols = cols]
  
  ## get every 3rd row of data talble
  codon_dt <- grdf[(seq(nrow(grdf)) - 1) %% 3 == 0, 1:5]
  #codon_dt$codon_number <-  unlist(sapply(txlen$cds_len, function(x){1:(x/3)}))
  
  #ribo_nt <- grdf[, .(riboseq = sum(riboseq), codon = paste(nt, sep="", collapse="")), by= (seq(nrow(grdf)) - 1) %/% 3]
  cols <- c("riboseq", paste0("X", riboseq_lengths))
  ribo_len <- grdf[, lapply(.SD, sum), by = (seq(nrow(grdf)) - 1) %/% 3, .SDcols = cols]
  codon_dt <- cbind(codon_dt, ribo_len)
  #codons <- grdf[, .(codon = paste(nt, sep="", collapse="")), by= (seq(nrow(grdf)) - 1) %/% 3]
  
  ## well expressed transcripts (cut-off) / with big distribution of footprints (median on codon coverage)
  rssum <- codon_dt[, .(sum = sum(riboseq)), by = "seqnames"]
  #rsmedian <- codon_dt[, .(median = median(riboseq)), by = "seqnames"]
  wet <- rssum[rssum$sum > 300,]$seqnames
  
  ## subset wet, convert into zscores
  #wet_gr <- gr[seqnames(gr) %in% wet,]
  wet_codon_dt <- codon_dt[seqnames %in% wet,]
  
  ## do not use 1st and last codon for zscore calculation
  #wet_txlen <- txlen[rownames(txlen) %in% wet,]
  #ends <- cumsum(wet_txlen$cds_len / 3)
  #starts <- c(1, ends[1:(length(ends)-1)]+1)
  #wet_codon_dt_inside <- wet_codon_dt[-c(starts, ends),]
  
  ## convert into zscores
  #zs <- wet_codon_dt_inside[, .(zscore = as.numeric(scale(riboseq))), by = "seqnames"]
  #wet_codon_dt_inside$zscore <- zs$zscore
  zs <- wet_codon_dt[, .(zscore = as.numeric(scale(riboseq))), by = "seqnames"]
  wet_codon_dt$zscore <- zs$zscore
  
  ## call peaks (stall sites)
  #peaks <- wet_codon_dt_inside[zscore > 8,]
  peaks <- wet_codon_dt[zscore > 8,]
  
  ### save
  save(peaks, file = c(file.path(peaks_path, libs[i])))
}


### load peaks for each study, get consensus (find overlaps)
peaks_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/peaks"
consensus_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/consensus"
libs <- list.files(path = peaks_path)
studies <- unique(sapply(libs, function(x){strsplit(x, "_")[[1]][1]}))

for (i in 1:length(studies)) {
  pattern <- paste0(studies[i], ".+")
  runs <- list.files(path = peaks_path, pattern = pattern)
  dt <- data.table(seqnames=c("transcript"), start=c(0))
  for (j in 1:length(runs)){
    load(file = c(file.path(peaks_path, runs[j]))) # peaks
    peaks <- peaks[,1:2]
    dt <- rbind(dt, peaks)
  }
  consensus <- count(dt, vars = c("seqnames", "start")) ## from library(plyr)
  save(consensus, file = c(file.path(consensus_path, paste0(studies[i], ".Rsave"))))
}



### get peaks with min 2 replicates
## sequence logos
libs <- list.files(path = consensus_path)
logos_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/logos"
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")

for (i in 1:length(libs)) {
  load(file = c(file.path(consensus_path, libs[i])))
  ss <- consensus[consensus$freq > 1,] # except for subtelny
  seqs <- c()
  for (j in 1:nrow(ss)) {
    tx <- as.character(ss[j,]$seqnames)
    start <- ss[j,]$start
    if (start > 21) {
      seq <- toupper(paste(fasta_cdna[[tx]][(start-21):(start+23)], sep="", collapse=""))
      seqs <- c(seqs, seq)
    }
  }
  ## write seqs to text file
  fileConn <- file(c(file.path(logos_path, paste0(substr(libs[i], 1, nchar(libs[i])-6), ".txt"))))
  writeLines(seqs, fileConn)
  close(fileConn)
}


### plot features:
# table with footprint lengths around peaks (for footprint shift plot)

### 'conservation' in tissues: group by tissue, find overlaps between studies -> UpSetR

### compare to peaks in other organisms

### stall sites per codon - what are the most frequent? P-site & A-site

### 18S complementarity: 2016_03_14_labmeeting

#####################
### stall sites in uORFs
### stall sites in 3'trailers

### group by tissue
HEK293 <- c("andreev2015.Rsave", "ingolia2012.Rsave", "ingolia2014.Rsave", "lee2012.Rsave",
            "liu2013.Rsave", "sidrauski2015.Rsave", "subtelny2014.Rsave")
THP1 <- c("fritsch2012.Rsave")
brain <- c("gonzalez2014.Rsave")
HeLa <- c("guo2010.Rsave", "stumpf2013.Rsave", "yoon2014.Rsave") ## remove "yoon2014.Rsave"
PC3 <- c("hsieh2012.Rsave")
fibroblasts <- c("rutkowski2015.Rsave", "stern2012.Rsave")

####################

## peaks <- peaks[riboseq > 30,] ## ?????????????????????????????????????????????????????????

### consensus from all libs
peaks_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/peaks_excl_15_5"
libs <- list.files(path = peaks_path)

dt <- data.table(seqnames=c("transcript"), start=c(0), library=c("lib"))
for (i in 1:length(libs)) {
  load(file = c(file.path(peaks_path, libs[i]))) # peaks
  peaks <- peaks[riboseq > 30,] ### get peaks higher than 30
  peaks <- peaks[,1:2]
  peaks$library <- rep(substr(libs[i], 1, nchar(libs[i])-6), nrow(peaks))
  dt <- rbind(dt, peaks)
}
consensus <- count(dt, vars = c("seqnames", "start")) ## from library(plyr)
save(consensus, file = "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/consensus_all_excl_15_5.Rsave")
save(dt, file = "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/all_peaks_excl_15_5.Rsave")

ggplot(consensus, aes(x = freq)) + geom_histogram()


ss <- consensus[consensus$freq > 4,]
seqs <- c()
snames <- c()
for (j in 1:nrow(ss)) {
  tx <- as.character(ss[j,]$seqnames)
  start <- ss[j,]$start
  if (start > 21) {
    seq <- toupper(paste(fasta_cdna[[tx]][(start-21):(start+23)], sep="", collapse=""))
    seqs <- c(seqs, seq)
    sname <- paste0(tx, "_", start)
    snames <- c(snames, sname)
  }
}
names(seqs) <- snames
## write seqs to text file
snames <- snames[nchar(seqs) == 45]
seqs <- seqs[nchar(seqs) == 45]

######### without PROLINES CC*
seqs <- seqs[substr(seqs, 22, 23) != "CC"] # at P-site
seqs <- seqs[substr(seqs, 25, 26) != "CC"] # at A-site
seqs <- seqs[substr(seqs, 19, 20) != "CC"] # at E-site

seqs <- lapply(seqs, as.SeqFastadna)
#names(seqs) <- snames

fasta_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/fasta"
fasta_file <- file.path(fasta_path, "consensus_all_4.txt")
write.fasta(seqs, snames, file.out = fasta_file, open = "w", nbchar = 60, as.string = FALSE)

seqs <- sample(seqs, 1000)
fasta_file <- file.path(fasta_path, "consensus_all_10_excl_15_5.txt")
write.fasta(seqs, names(seqs), file.out = fasta_file, open = "w", nbchar = 60, as.string = FALSE)

seqs <- lapply(seqs, function(x){substr(x, 7, 39)})
fasta_file <- file.path(fasta_path, "consensus_all_excl_15_5_thin.txt")
write.fasta(seqs, names(seqs), file.out = fasta_file, open = "w", nbchar = 60, as.string = FALSE)

fileConn <- file(c(file.path(logos_path, "consensus_all_10.txt")))
writeLines(seqs, fileConn)
close(fileConn)

##################################################################################################################
##### control / background
bcg <- txlen[rownames(txlen) %in% as.character(unique(consensus$seqnames)),]
bcg <- apply(bcg, 1, function(x){seq(as.numeric(x[7])+21, (as.numeric(x[7])+as.numeric(x[6])-23))})

for (i in 1:length(bcg)) {
  tx <- names(bcg[i])
  peaks <- consensus[consensus$seqnames == tx,]$start
  peak_region <- c()
  for (j in 1:length(peaks)) {
    region <- seq(peaks[j]-44, peaks[j]+64)
    peak_region <- c(peak_region, region)
  }
  peak_region <- unique(peak_region)
  bcg[[i]] <- bcg[[i]][! bcg[[i]] %in% peak_region]
}

sample_bcg <- sapply(bcg, function(x){sample(x,1)})

for (i in 1:length(sample_bcg)) {
  tx <- names(sample_bcg[i])
  start <- as.numeric(sample_bcg[i])
  if (start > 21) {
    seq <- toupper(paste(fasta_cdna[[tx]][(start-21):(start+23)], sep="", collapse=""))
    seqs <- c(seqs, seq)
    sname <- paste0(tx, "_", start)
    snames <- c(snames, sname)
  }
}
## write seqs to text file
snames <- snames[nchar(seqs) == 45]
seqs <- seqs[nchar(seqs) == 45]
seqs <- lapply(seqs, as.SeqFastadna)
names(seqs) <- snames

fasta_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/fasta"
fasta_file <- file.path(fasta_path, "background_all.txt")
write.fasta(seqs, snames, file.out = fasta_file, open = "w", nbchar = 60, as.string = FALSE)

##
seqs320 <- sample(seqs, 320)
fasta_file320 <- file.path(fasta_path, "background_sample320.txt")
write.fasta(seqs320, snames, file.out = fasta_file320, open = "w", nbchar = 60, as.string = FALSE)

##################################################################################################################



######### for peak in most libs
tx <- as.character(consensus[which.max(consensus$freq),]$seqnames)
temp <- dt[dt$seqnames %in% ss$seqnames,]
temp <- temp[temp$start %in% ss$start,]
temp <- temp[order(temp$seqnames, temp$start)]

libs <- temp[temp$seqnames == tx,]$library
libs <- as.character(sapply(libs, function(x){paste(x, ".Rsave", sep="", collapse="")}))
granges_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges"

aroundSS <- GRanges(Rle(tx, rep(1, 1)), IRanges((136-21), width=rep(45, 1)))

for (i in 1:length(libs)) {
  load(file = c(file.path(granges_path, libs[i]))) # peaks
  grSS <- subsetByOverlaps(gr, aroundSS, ignore.strand=TRUE)
  save(grSS, file = c(file.path("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/ENST00000262982", libs[i])))
}


########### make metaplots of fragment lengths
path_gr_ENST00000262982 <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/ENST00000262982"
libs <- list.files(path = path_gr_ENST00000262982)
save_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/PLOTS/ENST00000262982"

for (i in 1:length(libs)) {
  load(file = c(file.path(path_gr_ENST00000262982, libs[i]))) # peaks
  dt <- as.data.frame(grSS)
  setDT(dt)
  riboseq_lengths <- colnames(dt)[7:length(colnames(dt))]
  
  colourCount <- length(riboseq_lengths) # number of levels
  getPalette <- colorRampPalette(brewer.pal(11, "RdGy"))
  
  d <- dt[,7:length(colnames(dt))]
  d$scale <- c(-21:23)
  d <- melt(d, id.var="scale")
  p <- ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdGy"))(colourCount)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(fill=guide_legend(title="footprint\nlength"))
  ggsave(file = c(file.path(save_path, paste0(substr(libs[i], 1, nchar(libs[i])-6), ".png"))), plot = p)
}


###########################################
###########################################
#dupa <- data.table(seqnames = c(rep("dupa", 5), rep("chuj", 5)), start = c(1:5, 1:5), end = c(1:5, 1:5), kurwa = c(1:10))
#sub <- data.table(seqnames = c("dupa", "chuj"), start = c(3, 2), end = c(5,4))
#setkey(sub, seqnames, start, end)
#matches <- foverlaps(dupa,sub,type="within",nomatch=0L)

###########################################
all_peaks <- get(load(file = "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/all_peaks.Rsave")) # dt
all_peaks <- all_peaks[order(seqnames, start),]

peak_count <- count(all_peaks, vars = c("seqnames", "start"))
peak_count <- data.table(seqnames = rep(peak_count$seqnames, peak_count$freq),
           start = rep(peak_count$start, peak_count$freq), freq = rep(peak_count$freq, peak_count$freq))
peak_count$library <- all_peaks$library

# get peaks that exist in more than 10 libraries
peak10 <- peak_count[peak_count$freq > 10,]

t <- table(peak10$library)
runs <- names(t[t > 100])

granges_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges"
save_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/PLOTS/fragment_length"

for (i in 1:length(runs)) {
  peaks <- peak10[peak10$library == runs[i],]
  peaks <- peaks[start > 21,]
  aroundSS <- GRanges(Rle(peaks$seqnames, rep(1, 1)), IRanges((peaks$start-21), width=rep(45, 1)))
  
  load(file = c(file.path(granges_path, paste0(runs[i], ".Rsave")))) # gr
  
  stall_sites <- subsetByOverlaps(gr, aroundSS, type = "any")
  stall_sites <- split(stall_sites, seqnames(stall_sites))
  stall_sites <- stall_sites[sapply(stall_sites, function(x){length(x) == 45})]
  
  ss <- as.data.table(stall_sites[[1]])
  for (j in 2:length(stall_sites)) {
    ss <- ss + as.data.table(stall_sites[[j]])
  }
  
  ss <- ss[,7:ncol(ss)]
  riboseq_lengths <- colnames(ss)
  
  colourCount <- length(riboseq_lengths) # number of levels
  getPalette <- colorRampPalette(brewer.pal(11, "RdGy"))
  
  d <- ss
  d$scale <- c(-21:23)
  d <- melt(d, id.var="scale")
  p <- ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdGy"))(colourCount)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(fill=guide_legend(title="footprint\nlength"))
  ggsave(file = c(file.path(save_path, paste0(runs[i], ".png"))), plot = p)
}

### control fragment length




#################################################################
######## check for prolines and stop codons



