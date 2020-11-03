### plot distribution of means (on all codons of all txs included in stalling analysis)
### get bottom and top threshold (distribution tails), say 27.28 & 29.12
### plot these as vertical lines on stall site mean plots

library(ggplot2)
library(reshape2)
library(scales)

library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(RColorBrewer)
library(colorRamps)

library(seqinr)


### STOP CODON ###

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))

gr_path <- "/Volumes/USELESS/STALLING/gr"
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
peaks_path <- "/Volumes/USELESS/STALLING/peaks/median" # median or threshold, pick
save_stop_path <- "/Volumes/USELESS/STALLING/PLOTS/stop_avg"
save_stall_path <- "/Volumes/USELESS/STALLING/PLOTS/stall_avg"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  
  txLengths <- txLengths[order(rownames(txLengths)),]
  txLengths$utr5_cds_len <- txLengths$utr5_len + txLengths$cds_len
  # get longest transcript per gene
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  txlen <- txlen[txlen$cds_len > 100,]
  txlen <- txlen[order(rownames(txlen)),]
  
  if (org == "yeast") {
    # whole tx
    txlen_stop <- txlen
    cds_coord <- mapply(function(x,y){(x+4):(y-6)}, txlen_stop$utr5_len, txlen_stop$utr5_cds_len)
    tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord_wholetx <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_wholetx <- as.numeric(unlist(grdf_coord_wholetx))
    # stop
    txlen_stop <- txlen
    cds_coord <- sapply(txlen_stop$utr5_cds_len, function(x){(x-29):x})
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord_stop <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_stop <- split(grdf_coord_stop, rep(1:ncol(grdf_coord_stop), each = nrow(grdf_coord_stop)))
    grdf_coord_stop <- as.numeric(unlist(grdf_coord_stop))
  } else {
    # whole tx
    txlen_stop <- txlen
    cds_coord <- mapply(function(x,y){(x+4):(y-6)}, txlen_stop$utr5_len, txlen_stop$utr5_cds_len)
    tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord_wholetx <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_wholetx <- as.numeric(unlist(grdf_coord_wholetx))
    # stop
    txlen_stop <- txlen[txlen$utr3_len > 30,]
    cds_coord <- sapply(txlen_stop$utr5_cds_len, function(x){(x-17):(x+15)})
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen_stop)-1)]))
    grdf_coord_stop <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_stop <- split(grdf_coord_stop, rep(1:ncol(grdf_coord_stop), each = nrow(grdf_coord_stop)))
    grdf_coord_stop <- as.numeric(unlist(grdf_coord_stop))
  }
  
  
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)){
    load(file = c(file.path(gr_path, libs[j]))) # gr
    DT <- as.data.frame(gr)
    setDT(DT)
    
    # whole tx
    DT_wholetx <- DT[seqnames %in% txlen_stop$tx_name,]
    DT_wholetx <- DT_wholetx[grdf_coord_wholetx,]
    
    # stop
    DT_stop <- DT[seqnames %in% txlen_stop$tx_name,]
    DT_stop <- DT_stop[grdf_coord_stop,]
    
    # stall
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- peaks[start > 15]
    # for each peak, extract -15 + 17
    DT_stall <- DT[seqnames == as.character(peaks[1]$seqnames)][(peaks[1]$start-15):(peaks[1]$start+17)]
    for (r in 2:nrow(peaks)) {
      chunk <- DT[seqnames == as.character(peaks[r]$seqnames)][(peaks[r]$start-15):(peaks[r]$start+17)]
      if (all(!is.na(chunk))) {
        DT_stall <- rbind(DT_stall, chunk)
      }
    }
    
    ### distribution of means of DT_wholetx (per codon)
    # add up every 3 nt
    # DT_wholetx$codon <- rep(c(1:(nrow(DT_wholetx)/3)), each=3)
    # fld_wholetx <- DT_wholetx[, lapply(.SD, sum), by=codon, .SDcols=colnames(DT_wholetx)[7:(ncol(DT_wholetx)-1)] ]
    # colnames(fld_wholetx)[2:ncol(fld_wholetx)] <- sapply(colnames(fld_wholetx)[2:ncol(fld_wholetx)], function(x){substr(x,2,3)})
    # # multiply by fragment length
    # fl <- as.numeric(colnames(fld_wholetx[,2:(ncol(fld_wholetx)-1)]))
    # fld_wholetx$sum <- rowSums(as.data.frame(fld_wholetx)[,2:ncol(fld_wholetx)])
    # dupa <- t(apply(fld_wholetx[,2:(ncol(fld_wholetx)-1)], 1, function(x){x * fl}))
    # dupa <- as.data.frame(dupa)
    # setDT(dupa)
    # dupa$sum <- rowSums(as.data.frame(dupa)[,1:ncol(dupa)])
    # dupa$nof <- fld_wholetx$sum
    # dupa$avg <- dupa$sum / dupa$nof
    
    ## per nt
    colnames(DT_wholetx)[7:ncol(DT_wholetx)] <- sapply(colnames(DT_wholetx)[7:ncol(DT_wholetx)], function(x){substr(x,2,3)})
    fl <- as.numeric(colnames(DT_wholetx[,7:ncol(DT_wholetx)]))
    DT_wholetx$sum <- rowSums(as.data.frame(DT_wholetx)[,7:ncol(DT_wholetx)])
    dupa <- t(apply(DT_wholetx[,7:(ncol(DT_wholetx)-1)], 1, function(x){x * fl}))
    dupa <- as.data.frame(dupa)
    setDT(dupa)
    dupa$sum <- rowSums(as.data.frame(dupa)[,1:ncol(dupa)])
    dupa$nof <- DT_wholetx$sum
    dupa$avg <- dupa$sum / dupa$nof
    
    # ggplot(dupa, aes(avg)) + geom_density()  # distribution of avg
    # pick codons with X reads
    dupa2 <- dupa[dupa$nof > 10,]
    ggplot(dupa2, aes(avg)) + geom_density()
  
    # get tails and keep for plotting vertical lines (from dupa2)
    lower <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
    upper <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
    lower1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
    upper1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
    
    # plot per footprint length (from -10 nt if utr5 is available)
    if (org == "yeast") {
      DT_stop$group <- -30:-1
    } else {
      DT_stop$group <- c(-18:15)[-19]
    }
    
    DT_stall$group <- c(-15:18)[-16]
    
    fld_stop <- DT_stop[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stop)[7:(ncol(DT_stop)-1)] ]
    colnames(fld_stop)[2:ncol(fld_stop)] <- sapply(colnames(fld_stop)[2:ncol(fld_stop)], function(x){substr(x,2,3)})
    n <- ncol(fld_stop)
    fld_stop$avg <- apply(fld_stop[,2:n], 1, function(x){ sum(as.numeric(colnames(fld_stop[,2:n])) * x) / sum(x) })
    
    fld_stall <- DT_stall[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stall)[7:(ncol(DT_stall)-1)] ]
    colnames(fld_stall)[2:ncol(fld_stall)] <- sapply(colnames(fld_stall)[2:ncol(fld_stall)], function(x){substr(x,2,3)})
    n <- ncol(fld_stall)
    fld_stall$avg <- apply(fld_stall[,2:n], 1, function(x){ sum(as.numeric(colnames(fld_stall[,2:n])) * x) / sum(x) })
    
    # if (org == "yeast") {
    #   fld_stop$codon <- rep(c(-10:-1), each=3)
    #   fld_stop <- fld_stop[, mean(avg), by = c("codon")]
    # } else {
    #   fld_stop$codon <- rep(c(-5:6)[-6], each=3)
    #   fld_stop <- fld_stop[, mean(avg), by = c("codon")]
    # }
    
    if (org == "yeast") {
      fld_stop$codon <- c(-30:-1)
      fld_stop <- fld_stop[, mean(avg), by = c("codon")]
    } else {
      fld_stop$codon <- c(-15:18)[-16]
      fld_stop <- fld_stop[, mean(avg), by = c("codon")]
    }
    
    fld_stall$codon <- c(-15:18)[-16]
    fld_stall <- fld_stall[, mean(avg), by = c("codon")]
    
    #fld_stop$codon <- rep(c(-5:6)[-6], each=3)
    #fld_stop <- fld_stop[, mean(avg), by = c("codon")]
    
    ggplot(fld_stop, aes(codon, V1)) + geom_line() +
      geom_hline(yintercept=lower, linetype="dashed", color = "red") +
      geom_hline(yintercept=upper, linetype="dashed", color = "red") +
      geom_hline(yintercept=lower1, linetype="dashed", color = "red") +
      geom_hline(yintercept=upper1, linetype="dashed", color = "red")
    
    ## save plot
    ggsave(file = c(file.path(save_stop_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
    ggplot(fld_stall, aes(codon, V1)) + geom_line() +
      geom_hline(yintercept=lower, linetype="dashed", color = "red") +
      geom_hline(yintercept=upper, linetype="dashed", color = "red") +
      geom_hline(yintercept=lower1, linetype="dashed", color = "red") +
      geom_hline(yintercept=upper1, linetype="dashed", color = "red")
    
    ## save plot
    ggsave(file = c(file.path(save_stall_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
  }
  
}


######### sum per codon
####### stall site



### plot distribution over averages for all nucleotides (codons?) in a transcript (exclude first and last)
### use tails to plot horizontal lines on the stop and stall plots

# fld_stop$codon <- rep(c(-5:6)[-6], each=3)
# dupa <- fld_stop[, mean(avg), by = c("codon")]
# ggplot(dupa, aes(codon, V1)) + geom_line()


