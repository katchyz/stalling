### find stall sites in yeast, fruit fly, zebrafish, mouse, human
library(seqinr)
library(GenomicFeatures)
library(dplyr)
library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(colorRamps)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
gr_path <- "/Volumes/USELESS/STALLING/gr"
#save_peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/threshold"
#save_peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/median"
save_peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
wet_path <- "/Volumes/USELESS/STALLING/peaks_all/wet"

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
           gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                   "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
           fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                     "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                     "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                     "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                     "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

#load(file.path(gr_path, "mouse_chx.Rsave"))

##### select organism #####
# org <- "mouse"
# org <- "yeast"
# org <- "fruitfly"
# org <- "zebrafish"
# org <- "human"
###########################
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

for (i in 1:length(organisms)) {
  org <- organisms[i]

  fasta_cds <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
  }
  
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  # get longest transcript per gene
  ltpg <- arrange(txLengths, gene_id, desc(tx_len))
  ltpg <- ltpg[ltpg$cds_len > 0,]
  ltpg <- ltpg[!duplicated(ltpg$gene_id),]
  rownames(ltpg) <- ltpg$tx_name
  txlen <- txLengths
  # exons, cds
  exons <- exonsBy(txdb, by="tx", use.names=TRUE)
  cds <- cdsBy(txdb, by="tx", use.names=TRUE)
  
  ## for subsetting CDS by txlen coordinates 
  txlen <- txlen[txlen$cds_len %%3 == 0,] # get the cds which length is multiple of 3
  txlen <- txlen[txlen$tx_name %in% names(fasta_cds),] # get the tx that are in fasta
  txlen <- txlen[order(txlen$tx_name),]
  fasta_cds <- fasta_cds[names(fasta_cds) %in% txlen$tx_name]
  fasta_cds <- fasta_cds[order(names(fasta_cds))]
  fasta_cds <- fasta_cds[sapply(fasta_cds, function(x){length(x)}) == txlen$cds_len] ## check if fasta and txlen are equal
  txlen <- txlen[txlen$tx_name %in% names(fasta_cds),]
  
  ### cds_coord without first 15 and last 5 codons (Ingolia et al., 2011)
  # cds_coord <- mapply(function(x,y){(x+46):(x+y-15)}, txlen$utr5_len, txlen$cds_len)
  ### cds_coord without first and last coding codon
  cds_coord <- mapply(function(x,y){(x+4):(x+y-6)}, txlen$utr5_len, txlen$cds_len)
  tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
  grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
  grdf_coord <- unlist(grdf_coord)
  
  ############### LOAD GRANGES
  
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  
  for (i in 1:length(libs)) {
    ## load file
    load(file = c(file.path(gr_path, libs[i])))
    riboseq_lengths <- colnames(mcols(gr))[2:length(colnames(mcols(gr)))]
    grdf <- as.data.frame(gr)
    setDT(grdf)
    
    rm(gr)
    
    grdf <- grdf[seqnames %in% txlen$tx_name,] # get the longest CDSs, multiple of 3 only
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
    
    cols <- c("riboseq", paste0("X", riboseq_lengths))
    ribo_len <- grdf[, lapply(.SD, sum), by = (seq(nrow(grdf)) - 1) %/% 3, .SDcols = cols]
    codon_dt <- cbind(codon_dt, ribo_len)
    
    ## well expressed transcripts (cut-off) / with big distribution of footprints (median on codon coverage)
    #rssum <- codon_dt[, .(sum = sum(riboseq)), by = "seqnames"]
    rsmedian <- codon_dt[, .(median = median(riboseq)), by = "seqnames"]
    #wet <- rssum[rssum$sum > 300,]$seqnames
    wet <- rsmedian[rsmedian$median > 0,]$seqnames
    
    
    print(libs[i])
    print("Number of well-expressed genes:")
    print(length(unique(as.character(sapply(as.character(wet), function(x){txLengths[txLengths$tx_name == x,]$gene_id})))))
    
    
    ## subset wet, convert into zscores
    wet_codon_dt <- codon_dt[seqnames %in% wet,]
    
    ## convert into zscores
    zs <- wet_codon_dt[, .(zscore = as.numeric(scale(riboseq))), by = "seqnames"]
    wet_codon_dt$zscore <- zs$zscore
    
    save(wet_codon_dt, file = c(file.path(wet_path, libs[i])))
    
    ## call peaks (stall sites)
    #peaks <- wet_codon_dt_inside[zscore > 8,]
    peaks <- wet_codon_dt[zscore > 5,]
    
    ### save
    # save(peaks, file = c(file.path(save_peaks_path, libs[i])))
  }
}




#####################################################################################################
##### get consensus peaks for each organism #####

peaks_path <- "/Volumes/USELESS/STALLING/peaks/median"
#peaks_path <- "/Volumes/USELESS/STALLING/peaks/threshold"
consensus_path <- "/Volumes/USELESS/STALLING/consensus/median"
#consensus_path <- "/Volumes/USELESS/STALLING/consensus/threshold"

# org <- "yeast"
# org <- "fruitfly"
# org <- "zebrafish"
# org <- "mouse"
# org <- "human"

libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries

dt <- data.table(seqnames=c("transcript"), start=c(0))
for (i in 1:length(libs)){
  load(file = c(file.path(peaks_path, libs[i]))) # peaks
  peaks <- peaks[,1:2]
  dt <- rbind(dt, peaks)
}
consensus <- count(dt, vars = c("seqnames", "start")) ## from library(plyr)
save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))

# dt <- data.table(seqnames=c("transcript"), start=c(0))
# for (i in 1:length(libs)){
#   load(file = c(file.path(peaks_path, libs[i]))) # peaks
#   peaks <- peaks[peaks$zscore > 8]
#   peaks <- peaks[,1:2]
#   dt <- rbind(dt, peaks)
# }
# consensus <- count(dt, vars = c("seqnames", "start")) ## from library(plyr)
# sum(consensus$freq == 2)

##### if the tx is the same, position is off by 3 (adjacent codon) - merge / add up

off_peak <- 0
tx <- consensus[2,]$seqnames
start <- consensus[2,]$start
freq <- consensus[2,]$freq
for (i in 3:nrow(consensus)) {
  tx_new <- consensus[i,]$seqnames
  start_new <- consensus[i,]$start
  freq_new <- consensus[i,]$freq
  if (tx_new == tx) {
    if ((start_new - start) == 3) {
      off_peak <- off_peak + freq_new
    }
  }
  tx <- tx_new
  start <- start_new
  freq <- freq_new
}



######################################### z-score threshold = 8 ############################################
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/median"
#peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/threshold"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/median/8"
#consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/threshold/8"

peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
  dt <- data.table(seqnames=c("transcript"), start=c(0))
  for (j in 1:length(libs)){
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    #peaks <- peaks[peaks$zscore > 8]
    peaks <- peaks[,1:2]
    dt <- rbind(dt, peaks)
  }
  consensus <- plyr::count(dt, vars = c("seqnames", "start")) ## from library(plyr)
  save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))
}

########
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median_all_peaks"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
  dt <- data.table(seqnames=c("transcript"), start=c(0))
  for (j in 1:length(libs)){
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    #peaks <- peaks[peaks$zscore > 8]
    peaks <- peaks[,1:2]
    dt <- rbind(dt, peaks)
  }
  dt <- dt[2:nrow(dt),]
  consensus <- plyr::count(dt, vars = c("seqnames", "start")) ## from library(plyr)
  save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))
}


####################################################################################################
###### consensus - one transcript per gene #########################################################
####################################################################################################

###########
# load peaks, get only longest transcript, redo consensus
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/NEW_median"
organisms <- c("yeast", "fruitfly", "mouse", "human")

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
  dt <- data.table(seqnames=c("transcript"), start=c(0))
  ### longest tx
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  ##
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  ###
  for (j in 1:length(libs)){
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
    peaks <- peaks[,1:2]
    dt <- rbind(dt, peaks)
  }
  consensus <- plyr::count(dt, vars = c("seqnames", "start")) ## from library(plyr)
  setDT(consensus)
  consensus <- consensus[consensus$start > 15]
  consensus <- consensus[consensus$freq > 1]
  save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))
  #
  print(org)
  print(sum(consensus$freq == length(libs)))
}

# zebrafish:
# zebrafish_2h
# zebrafish_beaudoin_2h
# zebrafish_dome, zebrafish_shield
# zebrafish_subtelny_2h, zebrafish_subtelny_4h, zebrafish_subtelny_6h

org <- "zebrafish"
libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
dt <- data.table(seqnames=c("transcript"), start=c(0))
### longest tx
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# strand
strand <- transcripts(txdb)
strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
##
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
###
for (j in 1:2){
  load(file = c(file.path(peaks_path, libs[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt <- rbind(dt, peaks)
}
##
dt_chew <- data.table(seqnames=c("transcript"), start=c(0))
for (j in 3:4){
  load(file = c(file.path(peaks_path, libs[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt_chew <- rbind(dt_chew, peaks)
}
consensus_chew <- plyr::count(dt_chew, vars = c("seqnames", "start"))
consensus_chew <- consensus_chew[consensus_chew$freq > 1,]
setDT(consensus_chew)
consensus_chew <- consensus_chew[,1:2]
dt <- rbind(dt, consensus_chew)
##
dt_subtelny <- data.table(seqnames=c("transcript"), start=c(0))
for (j in 5:7){
  load(file = c(file.path(peaks_path, libs[j]))) # peaks
  peaks <- peaks[peaks$seqnames %in% rownames(txlen)] # one tx per gene
  peaks <- peaks[,1:2]
  dt_subtelny <- rbind(dt_subtelny, peaks)
}
consensus_subtelny <- plyr::count(dt_subtelny, vars = c("seqnames", "start"))
consensus_subtelny <- consensus_subtelny[consensus_subtelny$freq > 1,]
setDT(consensus_subtelny)
consensus_subtelny <- consensus_subtelny[,1:2]
dt <- rbind(dt, consensus_subtelny)
######
consensus <- plyr::count(dt, vars = c("seqnames", "start")) ## from library(plyr)
setDT(consensus)
consensus <- consensus[consensus$start > 15]
consensus <- consensus[consensus$freq > 1]
save(consensus, file = c(file.path(consensus_path, paste0(org, ".Rsave"))))

print(org)
print(sum(consensus$freq == 4))




#####################################################################
############# plot number of peaks ##################################

library(UpSetR)
library(GenomicFeatures)
library(plyr)

# load all peak libs
# put in df for each organism

peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
## yeast
y1 <- get(load(c(file.path(peaks_path, "yeast_gb15_meio2.Rsave"))))
y2 <- get(load(c(file.path(peaks_path, "yeast_chx.Rsave"))))
y3 <- get(load(c(file.path(peaks_path, "yeast_none.Rsave"))))

ss1 <- paste0(y1$seqnames, "_", y1$start)
ss2 <- paste0(y2$seqnames, "_", y2$start)
ss3 <- paste0(y3$seqnames, "_", y3$start)

# y1
y <- data.frame(ss = ss1, Y1 = rep(1,length(ss1)), Y2 = rep(0,length(ss1)), Y3 = rep(0,length(ss1)))
# y2
y$Y2[which(y$ss %in% ss2)] <- 1
y <- rbind(y, data.frame(ss = ss2[!ss2 %in% y$ss],
                         Y1 = rep(0, length(ss2[!ss2 %in% y$ss])),
                         Y2 = rep(1, length(ss2[!ss2 %in% y$ss])),
                         Y3 = rep(0, length(ss2[!ss2 %in% y$ss]))))
# y3
y$Y3[which(y$ss %in% ss3)] <- 1
y <- rbind(y, data.frame(ss = ss3[!ss3 %in% y$ss],
                         Y1 = rep(0, length(ss3[!ss3 %in% y$ss])),
                         Y2 = rep(0, length(ss3[!ss3 %in% y$ss])),
                         Y3 = rep(1, length(ss3[!ss3 %in% y$ss]))))

# remove first 5 codons
y <- y[as.numeric(sapply(as.character(y$ss), function(x){strsplit(x, split = "_")[[1]][2]})) > 15,]

upset(y, order.by = "freq", text.scale = 2, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("Y1", "Y2", "Y3"), color = "orange", active = T),
                     list(query = intersects, params = list("Y1", "Y2"), color = "orange", active = T),
                     list(query = intersects, params = list("Y1", "Y3"), color = "orange", active = T),
                     list(query = intersects, params = list("Y2", "Y3"), color = "orange", active = T)))


## fruit fly
f1 <- get(load(c(file.path(peaks_path, "fruitfly_S2cell_150_B.Rsave"))))
f2 <- get(load(c(file.path(peaks_path, "fruitfly_luo2018_S2_DMSO.Rsave"))))
f3 <- get(load(c(file.path(peaks_path, "fruitfly_SRR3031135.Rsave"))))

ss1 <- paste0(f1$seqnames, "_", f1$start)
ss2 <- paste0(f2$seqnames, "_", f2$start)
ss3 <- paste0(f3$seqnames, "_", f3$start)

# f1
f <- data.frame(ss = ss1, F1 = rep(1,length(ss1)), F2 = rep(0,length(ss1)), F3 = rep(0,length(ss1)))
# f2
f$F2[which(f$ss %in% ss2)] <- 1
f <- rbind(f, data.frame(ss = ss2[!ss2 %in% f$ss],
                         F1 = rep(0, length(ss2[!ss2 %in% f$ss])),
                         F2 = rep(1, length(ss2[!ss2 %in% f$ss])),
                         F3 = rep(0, length(ss2[!ss2 %in% f$ss]))))
# f3
f$F3[which(f$ss %in% ss3)] <- 1
f <- rbind(f, data.frame(ss = ss3[!ss3 %in% f$ss],
                         F1 = rep(0, length(ss3[!ss3 %in% f$ss])),
                         F2 = rep(0, length(ss3[!ss3 %in% f$ss])),
                         F3 = rep(1, length(ss3[!ss3 %in% f$ss]))))

# remove first 5 codons
f$tx <- sapply(as.character(f$ss), function(x){strsplit(x, split = "_")[[1]][1]})
f$pos <- as.numeric(sapply(as.character(f$ss), function(x){strsplit(x, split = "_")[[1]][2]}))
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Drosophila_melanogaster.BDGP6.79.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

f <- f[f$tx %in% rownames(txlen),]
f$utr5_len <- sapply(f$tx, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
f$stall <- f$pos - f$utr5_len
f <- f[f$stall > 15,]
f[,5:8] <- NULL

upset(f, order.by = "freq", text.scale = 2, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("F1", "F2", "F3"), color = "orange", active = T),
                     list(query = intersects, params = list("F1", "F2"), color = "orange", active = T),
                     list(query = intersects, params = list("F1", "F3"), color = "orange", active = T),
                     list(query = intersects, params = list("F2", "F3"), color = "orange", active = T)))


## zebrafish
z1 <- get(load(c(file.path(peaks_path, "zebrafish_bazzini2014.Rsave"))))
z2 <- get(load(c(file.path(peaks_path, "zebrafish_beaudoin_2h.Rsave"))))
z3 <- get(load(c(file.path(peaks_path, "zebrafish_subtelny_2h.Rsave"))))
z4 <- get(load(c(file.path(peaks_path, "zebrafish_dome.Rsave"))))
z5 <- get(load(c(file.path(peaks_path, "zebrafish_subtelny_4h.Rsave"))))
z6 <- get(load(c(file.path(peaks_path, "zebrafish_shield.Rsave"))))
z7 <- get(load(c(file.path(peaks_path, "zebrafish_subtelny_6h.Rsave"))))

ss1 <- paste0(z1$seqnames, "_", z1$start)
ss2 <- paste0(z2$seqnames, "_", z2$start)
ss3 <- paste0(z3$seqnames, "_", z3$start)
ss4 <- paste0(z4$seqnames, "_", z4$start)
ss5 <- paste0(z5$seqnames, "_", z5$start)
ss6 <- paste0(z6$seqnames, "_", z6$start)
ss7 <- paste0(z7$seqnames, "_", z7$start)

# z1
z <- data.frame(ss = ss1, Z1 = rep(1,length(ss1)), Z2 = rep(0,length(ss1)), Z3 = rep(0,length(ss1)),
                Z4 = rep(0,length(ss1)), Z5 = rep(0,length(ss1)), Z6 = rep(0,length(ss1)), Z7 = rep(0,length(ss1)))
# z2
z$Z2[which(z$ss %in% ss2)] <- 1
z <- rbind(z, data.frame(ss = ss2[!ss2 %in% z$ss],
                         Z1 = rep(0, length(ss2[!ss2 %in% z$ss])),
                         Z2 = rep(1, length(ss2[!ss2 %in% z$ss])),
                         Z3 = rep(0, length(ss2[!ss2 %in% z$ss])),
                         Z4 = rep(0, length(ss2[!ss2 %in% z$ss])),
                         Z5 = rep(0, length(ss2[!ss2 %in% z$ss])),
                         Z6 = rep(0, length(ss2[!ss2 %in% z$ss])),
                         Z7 = rep(0, length(ss2[!ss2 %in% z$ss]))))
# z3
z$Z3[which(z$ss %in% ss3)] <- 1
z <- rbind(z, data.frame(ss = ss3[!ss3 %in% z$ss],
                         Z1 = rep(0, length(ss3[!ss3 %in% z$ss])),
                         Z2 = rep(0, length(ss3[!ss3 %in% z$ss])),
                         Z3 = rep(1, length(ss3[!ss3 %in% z$ss])),
                         Z4 = rep(0, length(ss3[!ss3 %in% z$ss])),
                         Z5 = rep(0, length(ss3[!ss3 %in% z$ss])),
                         Z6 = rep(0, length(ss3[!ss3 %in% z$ss])),
                         Z7 = rep(0, length(ss3[!ss3 %in% z$ss]))))

# z4
z$Z4[which(z$ss %in% ss4)] <- 1
z <- rbind(z, data.frame(ss = ss4[!ss4 %in% z$ss],
                         Z1 = rep(0, length(ss4[!ss4 %in% z$ss])),
                         Z2 = rep(0, length(ss4[!ss4 %in% z$ss])),
                         Z3 = rep(0, length(ss4[!ss4 %in% z$ss])),
                         Z4 = rep(1, length(ss4[!ss4 %in% z$ss])),
                         Z5 = rep(0, length(ss4[!ss4 %in% z$ss])),
                         Z6 = rep(0, length(ss4[!ss4 %in% z$ss])),
                         Z7 = rep(0, length(ss4[!ss4 %in% z$ss]))))

# z5
z$Z5[which(z$ss %in% ss5)] <- 1
z <- rbind(z, data.frame(ss = ss5[!ss5 %in% z$ss],
                         Z1 = rep(0, length(ss5[!ss5 %in% z$ss])),
                         Z2 = rep(0, length(ss5[!ss5 %in% z$ss])),
                         Z3 = rep(0, length(ss5[!ss5 %in% z$ss])),
                         Z4 = rep(0, length(ss5[!ss5 %in% z$ss])),
                         Z5 = rep(1, length(ss5[!ss5 %in% z$ss])),
                         Z6 = rep(0, length(ss5[!ss5 %in% z$ss])),
                         Z7 = rep(0, length(ss5[!ss5 %in% z$ss]))))

# z6
z$Z6[which(z$ss %in% ss6)] <- 1
z <- rbind(z, data.frame(ss = ss6[!ss6 %in% z$ss],
                         Z1 = rep(0, length(ss6[!ss6 %in% z$ss])),
                         Z2 = rep(0, length(ss6[!ss6 %in% z$ss])),
                         Z3 = rep(0, length(ss6[!ss6 %in% z$ss])),
                         Z4 = rep(0, length(ss6[!ss6 %in% z$ss])),
                         Z5 = rep(0, length(ss6[!ss6 %in% z$ss])),
                         Z6 = rep(1, length(ss6[!ss6 %in% z$ss])),
                         Z7 = rep(0, length(ss6[!ss6 %in% z$ss]))))

# z7
z$Z7[which(z$ss %in% ss7)] <- 1
z <- rbind(z, data.frame(ss = ss7[!ss7 %in% z$ss],
                         Z1 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z2 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z3 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z4 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z5 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z6 = rep(0, length(ss7[!ss7 %in% z$ss])),
                         Z7 = rep(1, length(ss7[!ss7 %in% z$ss]))))


# remove first 5 codons
z$tx <- sapply(as.character(z$ss), function(x){strsplit(x, split = "_")[[1]][1]})
z$pos <- as.numeric(sapply(as.character(z$ss), function(x){strsplit(x, split = "_")[[1]][2]}))
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Danio_rerio.GRCz10.81_chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

z <- z[z$tx %in% rownames(txlen),]
z$utr5_len <- sapply(z$tx, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
z$stall <- z$pos - z$utr5_len
z <- z[z$stall > 15,]
z[,9:12] <- NULL

upset(z, sets = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"), order.by = "freq", text.scale = 1.5, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("Z1", "Z2"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z4", "Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z6", "Z4"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z5", "Z7", "Z3"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z2", "Z6", "Z4"), color = "orange", active = T),
                     list(query = intersects, params = list("Z5", "Z7", "Z6", "Z4"), color = "orange", active = T),
                     list(query = intersects, params = list("Z7", "Z3", "Z6", "Z4"), color = "orange", active = T)))


zb <- z[z$Z1 == 1,]
zb <- rbind(zb, z[z$Z2 == 1,])
zb <- rbind(zb, z[z$Z3 + z$Z5 + z$Z7 > 1,])
zb <- rbind(zb, z[z$Z4 + z$Z6 == 2,])

zbf <- zb[zb$Z1 + zb$Z2 + zb$Z3 + zb$Z4 + zb$Z5 + zb$Z6 + zb$Z7 > 1,]

# full
z$Z4_Z6 <- rep(0, nrow(z))
z$Z4_Z6[which(z$Z4 + z$Z6 == 2)] <- 1
z$Z3_Z5 <- rep(0, nrow(z))
z$Z3_Z5[which(z$Z3 + z$Z5 == 2)] <- 1
z$Z3_Z7 <- rep(0, nrow(z))
z$Z3_Z7[which(z$Z3 + z$Z7 == 2)] <- 1
z$Z5_Z7 <- rep(0, nrow(z))
z$Z5_Z7[which(z$Z5 + z$Z7 == 2)] <- 1
z$Z3 <- NULL
z$Z4 <- NULL
z$Z5 <- NULL
z$Z6 <- NULL
z$Z7 <- NULL

upset(z, sets = c("Z1", "Z2", "Z4_Z6", "Z3_Z5", "Z3_Z7", "Z5_Z7"), order.by = "freq", text.scale = 1.5, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("Z1", "Z2"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z3_Z5", "Z3_Z7", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z3_Z7", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z5_Z7", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z2", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z4_Z6", "Z3_Z5", "Z3_Z7", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z4_Z6", "Z3_Z5"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z3_Z5"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z3_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z2", "Z3_Z5"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z3_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z3_Z5", "Z3_Z7", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z2", "Z3_Z5", "Z3_Z7", "Z5_Z7"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z3_Z5", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z2", "Z3_Z5", "Z3_Z7", "Z5_Z7", "Z4_Z6"), color = "orange", active = T),
                     list(query = intersects, params = list("Z1", "Z3_Z5", "Z3_Z7", "Z5_Z7", "Z4_Z6"), color = "orange", active = T)))


## mouse
m1 <- get(load(c(file.path(peaks_path, "mouse_chx.Rsave"))))
m2 <- get(load(c(file.path(peaks_path, "mouse_none.Rsave"))))
m3 <- get(load(c(file.path(peaks_path, "mouse_3T3.Rsave"))))

ss1 <- paste0(m1$seqnames, "_", m1$start)
ss2 <- paste0(m2$seqnames, "_", m2$start)
ss3 <- paste0(m3$seqnames, "_", m3$start)

# m1
m <- data.frame(ss = ss1, M1 = rep(1,length(ss1)), M2 = rep(0,length(ss1)), M3 = rep(0,length(ss1)))
# m2
m$M2[which(m$ss %in% ss2)] <- 1
m <- rbind(m, data.frame(ss = ss2[!ss2 %in% m$ss],
                         M1 = rep(0, length(ss2[!ss2 %in% m$ss])),
                         M2 = rep(1, length(ss2[!ss2 %in% m$ss])),
                         M3 = rep(0, length(ss2[!ss2 %in% m$ss]))))
# m3
m$M3[which(m$ss %in% ss3)] <- 1
m <- rbind(m, data.frame(ss = ss3[!ss3 %in% m$ss],
                         M1 = rep(0, length(ss3[!ss3 %in% m$ss])),
                         M2 = rep(0, length(ss3[!ss3 %in% m$ss])),
                         M3 = rep(1, length(ss3[!ss3 %in% m$ss]))))

# remove first 5 codons
m$tx <- sapply(as.character(m$ss), function(x){strsplit(x, split = "_")[[1]][1]})
m$pos <- as.numeric(sapply(as.character(m$ss), function(x){strsplit(x, split = "_")[[1]][2]}))
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Mus_musculus.GRCm38.79.chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

m <- m[m$tx %in% rownames(txlen),]
m$utr5_len <- sapply(m$tx, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
m$stall <- m$pos - m$utr5_len
m <- m[m$stall > 15,]
m[,5:8] <- NULL

upset(m, order.by = "freq", text.scale = 2, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("M1", "M2", "M3"), color = "orange", active = T),
                     list(query = intersects, params = list("M1", "M2"), color = "orange", active = T),
                     list(query = intersects, params = list("M1", "M3"), color = "orange", active = T),
                     list(query = intersects, params = list("M2", "M3"), color = "orange", active = T)))


## human
h1 <- get(load(c(file.path(peaks_path, "human_fibroblasts_stern2012.Rsave"))))
h2 <- get(load(c(file.path(peaks_path, "human_stern_NODRUG.Rsave"))))
h3 <- get(load(c(file.path(peaks_path, "human_HeLa_stumpf2013.Rsave"))))
h4 <- get(load(c(file.path(peaks_path, "human_HEK_subtelny2014.Rsave"))))

ss1 <- paste0(h1$seqnames, "_", h1$start)
ss2 <- paste0(h2$seqnames, "_", h2$start)
ss3 <- paste0(h3$seqnames, "_", h3$start)
ss4 <- paste0(h4$seqnames, "_", h4$start)

# h1
h <- data.frame(ss = ss1, H1 = rep(1,length(ss1)), H2 = rep(0,length(ss1)), H3 = rep(0,length(ss1)),
                H4 = rep(0,length(ss1)))
# h2
h$H2[which(h$ss %in% ss2)] <- 1
h <- rbind(h, data.frame(ss = ss2[!ss2 %in% h$ss],
                         H1 = rep(0, length(ss2[!ss2 %in% h$ss])),
                         H2 = rep(1, length(ss2[!ss2 %in% h$ss])),
                         H3 = rep(0, length(ss2[!ss2 %in% h$ss])),
                         H4 = rep(0, length(ss2[!ss2 %in% h$ss]))))
# h3
h$H3[which(h$ss %in% ss3)] <- 1
h <- rbind(h, data.frame(ss = ss3[!ss3 %in% h$ss],
                         H1 = rep(0, length(ss3[!ss3 %in% h$ss])),
                         H2 = rep(0, length(ss3[!ss3 %in% h$ss])),
                         H3 = rep(1, length(ss3[!ss3 %in% h$ss])),
                         H4 = rep(0, length(ss3[!ss3 %in% h$ss]))))

# h4
h$H4[which(h$ss %in% ss4)] <- 1
h <- rbind(h, data.frame(ss = ss4[!ss4 %in% h$ss],
                         H1 = rep(0, length(ss4[!ss4 %in% h$ss])),
                         H2 = rep(0, length(ss4[!ss4 %in% h$ss])),
                         H3 = rep(0, length(ss4[!ss4 %in% h$ss])),
                         H4 = rep(1, length(ss4[!ss4 %in% h$ss]))))

# remove first 5 codons
h$tx <- sapply(as.character(h$ss), function(x){strsplit(x, split = "_")[[1]][1]})
h$pos <- as.numeric(sapply(as.character(h$ss), function(x){strsplit(x, split = "_")[[1]][2]}))
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

h <- h[h$tx %in% rownames(txlen),]
h$utr5_len <- sapply(h$tx, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
h$stall <- h$pos - h$utr5_len
h <- h[h$stall > 15,]
h[,6:9] <- NULL

upset(h, order.by = "freq", text.scale = 1.8, point.size = 5,
      mainbar.y.label = "No. peaks in common", sets.x.label = "No. peaks",
      queries = list(list(query = intersects, params = list("H1", "H2", "H3"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H2"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H3"), color = "orange", active = T),
                     list(query = intersects, params = list("H2", "H3"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H4"), color = "orange", active = T),
                     list(query = intersects, params = list("H2", "H4"), color = "orange", active = T),
                     list(query = intersects, params = list("H3", "H4"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H2", "H4"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H3", "H4"), color = "orange", active = T),
                     list(query = intersects, params = list("H1", "H2", "H3", "H4"), color = "orange", active = T)))


