### STALLING - analyze human libs
library(seqinr)
library(GenomicFeatures)
library(dplyr)
library(data.table)

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

cds_coord <- mapply(function(x,y){(x+1):(x+y)}, txlen$utr5_len, txlen$cds_len)
tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
grdf_coord <- unlist(grdf_coord)

#######################################################################################################
libs_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges"
peaks_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/peaks"
libs <- list.files(path = libs_path)
libs <- libs[64:96]

for (i in 1:length(libs)) {
  ## load file
  load(file = c(file.path(libs_path, libs[i])))
  #load("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges/andreev2015_SRR1173905.Rsave") # gr
  riboseq_lengths <- colnames(mcols(gr))[2:length(colnames(mcols(gr)))]
  grdf <- as.data.frame(gr)
  setDT(grdf)
  
  grdf <- grdf[seqnames %in% rownames(txlen),] # get the longest CDSs, multiple of 3 only
  grdf <- grdf[grdf_coord,] # get CDS coordinates
  
  ## add fasta... 
  grdf$nt <- as.character(unlist(fasta_cds))
  
  ## aggregate by codons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #df <- data.frame(seq=c('a','a','a','a','a','a','b','b','b'), riboseq = c(1,2,3,4,5,6,4,3,2), dupa = c(2,3,4,6,7,8,1,2,3), nt = c('a', 't', 'g', 'c', 'a', 'c', 'a', 't', 'g'))
  #setDT(df)
  #df[, .(riboseq = sum(riboseq), dupa = sum(dupa), codon = paste(nt, sep="", collapse="")), by= (seq(nrow(df)) - 1) %/% 3]
  #cols <- paste0("X", riboseq_lengths)
  #df[, lapply(.SD, sum), by = (seq(nrow(df)) - 1) %/% 3, .SDcols = cols]
  
  ## get every 3rd row of data talble
  codon_dt <- grdf[(seq(nrow(grdf)) - 1) %% 3 == 0, 1:5]
  codon_dt$codon_number <-  unlist(sapply(txlen$cds_len, function(x){1:(x/3)}))
  
  #ribo_nt <- grdf[, .(riboseq = sum(riboseq), codon = paste(nt, sep="", collapse="")), by= (seq(nrow(grdf)) - 1) %/% 3]
  cols <- c("riboseq", paste0("X", riboseq_lengths))
  ribo_len <- grdf[, lapply(.SD, sum), by = (seq(nrow(grdf)) - 1) %/% 3, .SDcols = cols]
  codon_dt <- cbind(codon_dt, ribo_len)
  #codons <- grdf[, .(codon = paste(nt, sep="", collapse="")), by= (seq(nrow(grdf)) - 1) %/% 3]
  
  ## well expressed transcripts (cut-off) / with big distribution of footprints (median on codon coverage)
  rssum <- codon_dt[, .(sum = sum(riboseq)), by = "seqnames"]
  #rsmedian <- codon_dt[, .(median = median(riboseq)), by = "seqnames"]
  wet <- rssum[rssum$sum > 500,]$seqnames
  
  ## subset wet, convert into zscores
  #wet_gr <- gr[seqnames(gr) %in% wet,]
  wet_codon_dt <- codon_dt[seqnames %in% wet,]
  
  ## do not use 1st and last codon for zscore calculation
  wet_txlen <- txlen[rownames(txlen) %in% wet,]
  ends <- cumsum(wet_txlen$cds_len / 3)
  starts <- c(1, ends[1:(length(ends)-1)]+1)
  
  wet_codon_dt_inside <- wet_codon_dt[-c(starts, ends),]
  
  ## convert into zscores
  zs <- wet_codon_dt_inside[, .(zscore = as.numeric(scale(riboseq))), by = "seqnames"]
  wet_codon_dt_inside$zscore <- zs$zscore
  
  ## call peaks (stall sites)
  peaks <- wet_codon_dt_inside[zscore > 8,]
  
  ### save
  save(peaks, file = c(file.path(peaks_path, libs[i])))
}


### load peaks for each study, get consensus (find overlaps)
### plot features:
# save nt sequence around peaks (for logo)
# table with footprint lengths around peaks (for footprint shift plot)
### 'conservation' in tissues: group by tissue, find overlaps between studies -> UpSetR

### compare to peaks in other organisms

### stall sites per codon - what are the most frequent? P-site & A-site

### 18S complementarity: 2016_03_14_labmeeting




