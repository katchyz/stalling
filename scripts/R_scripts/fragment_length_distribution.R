### fragment length distribution on CSSs (H1, H2, H3, H4)
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


H_table <- get(load("/Volumes/USELESS/STALLING/SNP/human_table_bam.Rsave"))
rm(h_tx_ss)

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]
txLengths_h$utr5_cds_len <- txLengths_h$utr5_len + txLengths_h$cds_len

gr_path <- "/Volumes/USELESS/STALLING/gr"
save_stop_path <- "/Volumes/USELESS/STALLING/PLOTS/stop_avg"
save_stall_path <- "/Volumes/USELESS/STALLING/PLOTS/stall_avg"
save_bckg_path <- "/Volumes/USELESS/STALLING/PLOTS/bckg_avg"

# whole tx
txlen_stop <- txLengths_h
cds_coord <- mapply(function(x,y){(x+4):(y-6)}, txlen_stop$utr5_len, txlen_stop$utr5_cds_len)
tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txLengths_h)-1)]))
grdf_coord_wholetx <- mapply(function(x,y){x + y}, cds_coord, tx_len)
grdf_coord_wholetx <- as.numeric(unlist(grdf_coord_wholetx))
# stop
txlen_stop <- txLengths_h[txLengths_h$utr3_len > 30,]
cds_coord <- sapply(txlen_stop$utr5_cds_len, function(x){(x-17):(x+15)})
cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
names(cds_coord) <- NULL
tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen_stop)-1)]))
grdf_coord_stop <- mapply(function(x,y){x + y}, cds_coord, tx_len)
grdf_coord_stop <- split(grdf_coord_stop, rep(1:ncol(grdf_coord_stop), each = nrow(grdf_coord_stop)))
grdf_coord_stop <- as.numeric(unlist(grdf_coord_stop))

org = "human"
libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
libs_order <- c("h4", "h3", "h1", "h2")

for (j in 1:length(libs)){
  load(file = c(file.path(gr_path, libs[j]))) # gr
  DT <- as.data.frame(gr)
  setDT(DT)
  rm(gr)
  
  # whole tx
  DT_wholetx <- DT[seqnames %in% txlen_stop$tx_name,]
  DT_wholetx <- DT_wholetx[grdf_coord_wholetx,]
  
  # stop
  DT_stop <- DT[seqnames %in% txlen_stop$tx_name,]
  DT_stop <- DT_stop[grdf_coord_stop,]
  
  # stall
  #load(file = c(file.path(peaks_path, libs[j]))) # peaks
  #peaks <- peaks[start > 15]
  lo <- libs_order[j]
  peaks <- H_table[H_table[lo] == 1,]
  # for each peak, extract -15 + 17
  #DT_stall <- DT[seqnames == as.character(peaks[1]$seqnames)][(peaks[1]$start-15):(peaks[1]$start+17)]
  DT_stall <- DT[seqnames == as.character(peaks[1,]$h_tx)][(peaks[1,]$h_ss-15):(peaks[1,]$h_ss+17)]
  for (r in 2:nrow(peaks)) {
    chunk <- DT[seqnames == as.character(peaks[r,]$h_tx)][(peaks[r,]$h_ss-15):(peaks[r,]$h_ss+17)]
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
  DT_wholetx <- DT_wholetx[DT_wholetx$riboseq > 0,]
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
  ggplot(dupa2, aes(avg)) + geom_density() + theme_classic() +
    geom_vline(xintercept=lower1, linetype="dashed", color = "red") +
    geom_vline(xintercept=upper1, linetype="dashed", color = "red")
  
  ggsave(file = c(file.path(save_bckg_path, paste0(libs_order[j], ".png"))), height = 5 , width = 7)
  
  # get tails and keep for plotting vertical lines (from dupa2)
  #lower <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
  #upper <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
  lower1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
  upper1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
  
  # plot per footprint length (from -10 nt if utr5 is available)
  DT_stop$group <- rep(c(-18:15)[-19], nrow(DT_stop)/33)
  DT_stall$group <- rep(c(-15:18)[-16], nrow(DT_stall)/33)
  
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
  

  fld_stop$codon <- c(-15:18)[-16]
  fld_stop <- fld_stop[, mean(avg), by = c("codon")]
  
  fld_stall$codon <- c(-15:18)[-16]
  fld_stall <- fld_stall[, mean(avg), by = c("codon")]
  
  #fld_stop$codon <- rep(c(-5:6)[-6], each=3)
  #fld_stop <- fld_stop[, mean(avg), by = c("codon")]
  
  ggplot(fld_stop, aes(codon, V1)) + geom_line() +
    geom_hline(yintercept=lower1, linetype="dashed", color = "red") +
    geom_hline(yintercept=upper1, linetype="dashed", color = "red") +
    theme_classic() + xlab("distance from stop [nt]") + ylab("avg fragment length")
  
  ## save plot
  ggsave(file = c(file.path(save_stop_path, paste0(libs_order[j], ".png"))), height = 5 , width = 7)
  
  ggplot(fld_stall, aes(codon, V1)) + geom_line() +
    geom_hline(yintercept=lower1, linetype="dashed", color = "red") +
    geom_hline(yintercept=upper1, linetype="dashed", color = "red") +
    theme_classic() + xlab("distance from stall [nt]") + ylab("avg fragment length")
  
  ## save plot
  ggsave(file = c(file.path(save_stall_path, paste0(libs_order[j], ".png"))), height = 5 , width = 7)
  
}

#####################################################
########### 13 NMD ##################################

save_stall_path <- "/Volumes/USELESS/STALLING/PLOTS/stall_avg/nmd13"

for (j in 1:length(libs)){
  load(file = c(file.path(gr_path, libs[j]))) # gr
  DT <- as.data.frame(gr)
  setDT(DT)
  rm(gr)
  
  # whole tx
  DT_wholetx <- DT[seqnames %in% txlen_stop$tx_name,]
  DT_wholetx <- DT_wholetx[grdf_coord_wholetx,]
  
  # stall
  lo <- libs_order[j]
  peaks <- nmd13
  # for each peak, extract -15 + 17
  #DT_stall <- DT[seqnames == as.character(peaks[1]$seqnames)][(peaks[1]$start-15):(peaks[1]$start+17)]
  DT_stall <- DT[seqnames == as.character(peaks[1,]$tx)][(peaks[1,]$ss-15):(peaks[1,]$ss+17)]
  for (r in 2:nrow(peaks)) {
    chunk <- DT[seqnames == as.character(peaks[r,]$tx)][(peaks[r,]$ss-15):(peaks[r,]$ss+17)]
    if (all(!is.na(chunk))) {
      DT_stall <- rbind(DT_stall, chunk)
    }
  }
  
  
  ## per nt
  colnames(DT_wholetx)[7:ncol(DT_wholetx)] <- sapply(colnames(DT_wholetx)[7:ncol(DT_wholetx)], function(x){substr(x,2,3)})
  fl <- as.numeric(colnames(DT_wholetx[,7:ncol(DT_wholetx)]))
  DT_wholetx <- DT_wholetx[DT_wholetx$riboseq > 0,]
  DT_wholetx$sum <- rowSums(as.data.frame(DT_wholetx)[,7:ncol(DT_wholetx)])
  dupa <- t(apply(DT_wholetx[,7:(ncol(DT_wholetx)-1)], 1, function(x){x * fl}))
  dupa <- as.data.frame(dupa)
  setDT(dupa)
  dupa$sum <- rowSums(as.data.frame(dupa)[,1:ncol(dupa)])
  dupa$nof <- DT_wholetx$sum
  dupa$avg <- dupa$sum / dupa$nof
  
  ggplot(dupa, aes(avg)) + geom_density()  # distribution of avg
  # pick codons with X reads
  dupa2 <- dupa[dupa$nof > 10,]
  # ggplot(dupa2, aes(avg)) + geom_density() + theme_classic() +
  #   geom_vline(xintercept=lower1, linetype="dashed", color = "red") +
  #   geom_vline(xintercept=upper1, linetype="dashed", color = "red")
  
  # ggsave(file = c(file.path(save_bckg_path, paste0(libs_order[j], ".png"))), height = 5 , width = 7)
  
  # get tails and keep for plotting vertical lines (from dupa2)
  #lower <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
  #upper <- qnorm(0.05, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
  lower1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = T)
  upper1 <- qnorm(0.1, mean = mean(dupa2$avg), sd = sd(dupa2$avg), lower.tail = F)
  
  # plot per footprint length (from -10 nt if utr5 is available)
  DT_stall$group <- rep(c(-15:18)[-16], nrow(DT_stall)/33)
  
  fld_stall <- DT_stall[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stall)[7:(ncol(DT_stall)-1)] ]
  colnames(fld_stall)[2:ncol(fld_stall)] <- sapply(colnames(fld_stall)[2:ncol(fld_stall)], function(x){substr(x,2,3)})
  n <- ncol(fld_stall)
  fld_stall$avg <- apply(fld_stall[,2:n], 1, function(x){ sum(as.numeric(colnames(fld_stall[,2:n])) * x) / sum(x) })
  
  fld_stall$codon <- c(-15:18)[-16]
  fld_stall <- fld_stall[, mean(avg), by = c("codon")]
  
  
  ggplot(fld_stall, aes(codon, V1)) + geom_line() +
    geom_hline(yintercept=lower1, linetype="dashed", color = "red") +
    geom_hline(yintercept=upper1, linetype="dashed", color = "red") +
    theme_classic() + xlab("distance from stall [nt]") + ylab("avg fragment length")
  
  ## save plot
  ggsave(file = c(file.path(save_stall_path, paste0(libs_order[j], ".png"))), height = 5 , width = 7)
  
}

# fld_stall_H4 <- fld_stall

fld_stall <- data.table(codon = c(-15:18)[-16],
                        V1 = (fld_stall_H1$V1 + fld_stall_H2$V1 + fld_stall_H3$V1 + fld_stall_H4$V1)/4)








