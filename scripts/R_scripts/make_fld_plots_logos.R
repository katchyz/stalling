### make fragment length distribution plots
### save logos
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(seqinr)
library(data.table)
library(GenomicFeatures)
library(plyr)
library(rtracklayer)

consensus_path <- "/Volumes/USELESS/STALLING/consensus/median/5"
gr_path <- "/Volumes/USELESS/STALLING/gr"
peaks_path <- "/Volumes/USELESS/STALLING/peaks_whole_tx"
save_fld_path <- "/Volumes/USELESS/STALLING/plots_whole_tx/fld"
save_logo_path <- "/Volumes/USELESS/STALLING/plots_whole_tx/logo"
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
save_logo_path_long <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall_long"
# org <- "yeast"
# org <- "fruitfly"
# org <- "zebrafish"
# org <- "mouse"
# org <- "human"

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

### fld
for (i in 1:length(organisms)) {
  org <- organisms[i]
  # load consensus peaks, get those present in at least 2 libraries
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus <- consensus[consensus$freq > 1,]
  setDT(consensus)
  consensus <- consensus[consensus$start > 15]
  
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  
  consensus$utr5_cds_len <- sapply(consensus$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len + txLengths[txLengths$tx_name == x,]$cds_len})
  consensus <- consensus[consensus$start < (consensus$utr5_cds_len-11)]
  
  # load granges for each library
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)){
    load(file = c(file.path(gr_path, libs[j]))) # gr
    DT <- as.data.frame(gr)
    setDT(DT)
    # check if the peaks are present in the given library (peaks Rsave file)
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- as.data.frame(peaks)
    setDT(peaks)
    consensus <- consensus[consensus[peaks, on=.(seqnames,start), which=TRUE, nomatch=0]]
    ## subset on 'seqnames' and 'start' columns from consensus, get 15nt around codon with peak
    peak_idx <- DT[consensus, on=.(seqnames,start), which=TRUE, nomatch=0]
    surr_idx <- lapply(peak_idx, function(x){(x-15):(x+17)})
    ssDT <- DT[unlist(surr_idx),]
    ## get every 33 rows, add up each fragment length
    ssDT$group <- 1:33
    fld <- ssDT[, lapply(.SD, sum, na.rm=TRUE), by=group, .SDcols=colnames(ssDT)[7:(ncol(ssDT)-1)] ]
    ## plot
    colnames(fld)[2:ncol(fld)] <- sapply(colnames(fld)[2:ncol(fld)], function(x){substr(x,2,3)})
    fld$scale <- c(-15:17)
    fld$group <- NULL
    d <- melt(fld, id.var="scale")
    ## add extra shades to the scale
    colourCount <- ncol(fld)-1 # number of levels
    getPalette <- colorRampPalette(brewer.pal(11, "RdGy"))
    ## plot
    ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdGy"))(colourCount)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength"))
    #ggplot(d, aes(x=scale, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + scale_fill_brewer(name=element_blank(), palette="RdGy") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(fill=guide_legend(title="footprint\nlength"))
    ## save plot
    ggsave(file = c(file.path(save_fld_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))))
  }
}

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))


save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/NEW_median"
### logos
for (i in 1:length(organisms)) {
  org <- organisms[i]
  # load consensus peaks, get those present in at least 2 libraries
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus <- consensus[consensus$freq > 1,]
  setDT(consensus)
  consensus <- consensus[consensus$start > 15]
  # get fasta file around these peaks
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-15):(y+17)]}, as.character(consensus$seqnames), consensus$start))
  logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
  logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
  # write to file
  fileConn<-file(c(file.path(save_logo_path, paste0(org, ".txt"))))
  writeLines(logo, fileConn)
  close(fileConn)
}


### logos - LONG
for (i in 1:length(organisms)) {
  org <- organisms[i]
  # load consensus peaks, get those present in at least 2 libraries
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  consensus <- consensus[consensus$freq > 1,]
  setDT(consensus)
  consensus <- consensus[consensus$start > 90]
  # get fasta file around these peaks
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-90):(y+17)]}, as.character(consensus$seqnames), consensus$start))
  logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
  logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
  # write to file
  fileConn<-file(c(file.path(save_logo_path_long, paste0(org, ".txt"))))
  writeLines(logo, fileConn)
  close(fileConn)
}


#####################################################################################################################
#####################################################################################################################
##### make 5'end plots to see how many codons to exclude from the start
save_5end_path <- "/Volumes/USELESS/STALLING/5end"
save_5end_fld_path <- "/Volumes/USELESS/STALLING/5end_fld"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  
  txLengths <- txLengths[order(rownames(txLengths)),]
  # get longest transcript per gene
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  txlen <- txlen[txlen$cds_len > 100,]
  txlen <- txlen[order(rownames(txlen)),]
  
  if (org == "yeast") {
    cds_coord <- sapply(txlen$utr5_len, function(x){(x+1):(x+100)}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord <- split(grdf_coord, rep(1:ncol(grdf_coord), each = nrow(grdf_coord)))
    grdf_coord <- as.numeric(unlist(grdf_coord))
  } else {
    txlen <- txlen[txlen$utr5_len > 10,]
    cds_coord <- sapply(txlen$utr5_len, function(x){(x-9):(x+100)}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord <- split(grdf_coord, rep(1:ncol(grdf_coord), each = nrow(grdf_coord)))
    grdf_coord <- as.numeric(unlist(grdf_coord))
  }
  
  
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)){
    load(file = c(file.path(gr_path, libs[j]))) # gr
    DT <- as.data.frame(gr)
    setDT(DT)
    
    DT <- DT[seqnames %in% txlen$tx_name,]
    DT <- DT[grdf_coord,] # get first 100nt of each CDS coordinates
    
    # get every 30 rows (or group by seqnames), sum up riboseq column
    # d <- DT[, lapply(riboseq, sum), by=seqnames ]
    # d$seqnames <- NULL
    # d <- colSums(d)
    # d <- data.frame(scale=1:100, reads=as.numeric(d))
    # 
    # ggplot(d, aes(scale,reads)) + geom_bar(stat="identity") + ggtitle(substr(libs[j],1,nchar(libs[j])-6))
    # ggsave(file = c(file.path(save_5end_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))))
    
    # plot per footprint length (from -10 nt if utr5 is available)
    if (org == "yeast") {
      DT$group <- 1:100
    } else {
      DT$group <- c(-10:100)[-11]
    }
    
    fld <- DT[, lapply(.SD, sum), by=group, .SDcols=colnames(DT)[7:(ncol(DT)-1)] ]
    colnames(fld)[2:ncol(fld)] <- sapply(colnames(fld)[2:ncol(fld)], function(x){substr(x,2,3)})
    d <- melt(fld, id.var="group")

    ## plot
    ggplot(d, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity") +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength"))
    ## save plot
    ggsave(file = c(file.path(save_5end_fld_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))))
    
  }
}

### RdYlBu

#####################################################################################################################
#####################################################################################################################
##### make 3'end plots
save_3end_fld_path <- "/Volumes/USELESS/STALLING/3end_fld"

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
    cds_coord <- sapply(txlen$utr5_cds_len, function(x){(x-99):x}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord <- split(grdf_coord, rep(1:ncol(grdf_coord), each = nrow(grdf_coord)))
    grdf_coord <- as.numeric(unlist(grdf_coord))
  } else {
    txlen <- txlen[txlen$utr3_len > 10,]
    cds_coord <- sapply(txlen$utr5_len, function(x){(x-99):(x+10)}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord <- split(grdf_coord, rep(1:ncol(grdf_coord), each = nrow(grdf_coord)))
    grdf_coord <- as.numeric(unlist(grdf_coord))
  }
  
  
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)){
    load(file = c(file.path(gr_path, libs[j]))) # gr
    DT <- as.data.frame(gr)
    setDT(DT)
    
    DT <- DT[seqnames %in% txlen$tx_name,]
    DT <- DT[grdf_coord,] # get first 100nt of each CDS coordinates
    
    # plot per footprint length (from -100 to +10 nt if utr3 is available)
    if (org == "yeast") {
      DT$group <- 1:100
    } else {
      DT$group <- c(-100:10)[-101]
    }
    
    fld <- DT[, lapply(.SD, sum), by=group, .SDcols=colnames(DT)[7:(ncol(DT)-1)] ]
    colnames(fld)[2:ncol(fld)] <- sapply(colnames(fld)[2:ncol(fld)], function(x){substr(x,2,3)})
    d <- melt(fld, id.var="group")
    
    ## plot
    ggplot(d, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity") +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength"))
    ## save plot
    ggsave(file = c(file.path(save_3end_fld_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))))
    
  }
}






