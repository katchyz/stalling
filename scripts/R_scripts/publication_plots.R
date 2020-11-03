### plots

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

### example ###

fld <- data.table(group=c(-2,-1,1,2,3,4,5), X21=c(sample(1:10,7,replace=T)), X22=sample(1:10,7,replace=T), X23=sample(1:10,7,replace=T),
                  X24=sample(1:10,7,replace=T), X25=sample(1:10,7,replace=T), X26=sample(1:10,7,replace=T),
                  X27=sample(1:10,7,replace=T), X28=sample(1:10,7,replace=T), X29=sample(1:10,7,replace=T),
                  X30=sample(1:10,7,replace=T), X31=sample(1:10,7,replace=T), X32=sample(1:10,7,replace=T),
                  X33=sample(1:10,7,replace=T), X34=sample(1:10,7,replace=T), X35=sample(1:10,7,replace=T))

colnames(fld)[2:ncol(fld)] <- sapply(colnames(fld)[2:ncol(fld)], function(x){substr(x,2,3)})
d <- melt(fld, id.var="group")

## add extra shades to the scale
colourCount <- ncol(fld)-1 # number of levels
getPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))
## plot
# ggplot(d, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
#  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
#  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
#  guides(fill=guide_legend(title="footprint\nlength")) +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"))


fld2 <- data.table(group=c(-5,-4,-3,-2,-1,1,2), X21=c(sample(1:20,7,replace=T)), X22=sample(1:20,7,replace=T), X23=sample(1:20,7,replace=T),
                  X24=sample(1:20,7,replace=T), X25=sample(1:20,7,replace=T), X26=sample(1:20,7,replace=T),
                  X27=sample(1:20,7,replace=T), X28=sample(1:20,7,replace=T), X29=sample(1:20,7,replace=T),
                  X30=sample(1:20,7,replace=T), X31=sample(1:20,7,replace=T), X32=sample(1:20,7,replace=T),
                  X33=sample(1:20,7,replace=T), X34=sample(1:20,7,replace=T), X35=sample(1:20,7,replace=T))

colnames(fld2)[2:ncol(fld2)] <- sapply(colnames(fld2)[2:ncol(fld2)], function(x){substr(x,2,3)})
d2 <- melt(fld2, id.var="group")

d$feature <- c(rep("5'UTR",2), rep("5' CDS",5))
d2$feature <- c(rep("3' CDS",5), rep("3'UTR",2))

D <- rbind(d,d2)
D$feature <- factor(D$feature, levels=c("5'UTR", "5' CDS", "3' CDS", "3'UTR"))

ggplot(D, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill=guide_legend(title="footprint\nlength")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(~ feature, scales = "free_x", space = "free_x") + scale_x_continuous(breaks= pretty_breaks(2)) +
  theme(panel.spacing.x=unit(c(1,3,1), "lines"))


ggsave("/Users/kasia/Desktop/dupa.png", height = 5 , width = 20)

#############################################################################################################
##### 5'UTR, 5'CDS, 3'CDS, 3'UTR (for yeast: CDS only) #####
### geom_bar(stat="identity")

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))

gr_path <- "/Volumes/USELESS/STALLING/gr"
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/cds_utr"

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
    # start
    txlen_start <- txlen
    cds_coord <- sapply(txlen_start$utr5_len, function(x){(x+1):(x+60)}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen_start$tx_len[1:(nrow(txlen_start)-1)]))
    grdf_coord_start <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_start <- split(grdf_coord_start, rep(1:ncol(grdf_coord_start), each = nrow(grdf_coord_start)))
    grdf_coord_start <- as.numeric(unlist(grdf_coord_start))
    # stop
    txlen_stop <- txlen
    cds_coord <- sapply(txlen_stop$utr5_cds_len, function(x){(x-59):x}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen_stop$tx_len[1:(nrow(txlen)-1)]))
    grdf_coord_stop <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_stop <- split(grdf_coord_stop, rep(1:ncol(grdf_coord_stop), each = nrow(grdf_coord_stop)))
    grdf_coord_stop <- as.numeric(unlist(grdf_coord_stop))
  } else {
    # start
    txlen_start <- txlen[txlen$utr5_len > 30,]
    cds_coord <- sapply(txlen_start$utr5_len, function(x){(x-29):(x+60)}) # rows are 1-100 coords
    cds_coord <- split(cds_coord, rep(1:ncol(cds_coord), each = nrow(cds_coord)))
    names(cds_coord) <- NULL
    tx_len <- cumsum(c(0, txlen_start$tx_len[1:(nrow(txlen_start)-1)]))
    grdf_coord_start <- mapply(function(x,y){x + y}, cds_coord, tx_len)
    grdf_coord_start <- split(grdf_coord_start, rep(1:ncol(grdf_coord_start), each = nrow(grdf_coord_start)))
    grdf_coord_start <- as.numeric(unlist(grdf_coord_start))
    # stop
    txlen_stop <- txlen[txlen$utr3_len > 30,]
    cds_coord <- sapply(txlen_stop$utr5_cds_len, function(x){(x-59):(x+30)}) # rows are 1-100 coords
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
    DT_start <- as.data.frame(gr)
    setDT(DT_start)
    DT_stop <- DT_start
    
    # start
    DT_start <- DT_start[seqnames %in% txlen_start$tx_name,]
    DT_start <- DT_start[grdf_coord_start,] 
    # stop
    DT_stop <- DT_stop[seqnames %in% txlen_stop$tx_name,]
    DT_stop <- DT_stop[grdf_coord_stop,]
    
    # plot per footprint length (from -10 nt if utr5 is available)
    if (org == "yeast") {
      DT_start$group <- 1:60
      DT_stop$group <- -60:-1
    } else {
      DT_start$group <- c(-30:60)[-31]
      DT_stop$group <- c(-60:30)[-61]
    }
    
    fld_start <- DT_start[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_start)[7:(ncol(DT_start)-1)] ]
    colnames(fld_start)[2:ncol(fld_start)] <- sapply(colnames(fld_start)[2:ncol(fld_start)], function(x){substr(x,2,3)})
    colourCount <- ncol(fld_start)-1 # add extra shades to the scale
    fld_start <- melt(fld_start, id.var="group")
    if (org == "yeast") {
      fld_start$feature <- rep("5' CDS",60)
    } else {
      fld_start$feature <- c(rep("5'UTR",30), rep("5' CDS",60))
    }
    
    
    fld_stop <- DT_stop[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stop)[7:(ncol(DT_stop)-1)] ]
    colnames(fld_stop)[2:ncol(fld_stop)] <- sapply(colnames(fld_stop)[2:ncol(fld_stop)], function(x){substr(x,2,3)})
    fld_stop <- melt(fld_stop, id.var="group")
    if (org == "yeast") {
      fld_stop$feature <- rep("3' CDS",60)
    } else {
      fld_stop$feature <- c(rep("3' CDS",60), rep("3'UTR",30))
    }
    
    # add extra shades to the scale
    # colourCount <- ncol(fld_start)-1 # number of levels
    getPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))
    
    D <- rbind(fld_start,fld_stop)
    if (org == "yeast") {
      D$feature <- factor(D$feature, levels=c("5' CDS", "3' CDS"))
      ## plot
      ggplot(D, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity") +
        scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background = element_blank(),
              strip.text.x = element_blank(), text = element_text(size=30),
              legend.justification = c(0.5, 0.9), legend.position = c(0.5, 0.9), legend.direction = "horizontal") +
        guides(fill=guide_legend(title="footprint\nlength")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        facet_grid(~ feature, scales = "free_x", space = "free_x") + scale_x_continuous(breaks= pretty_breaks(2)) +
        theme(panel.spacing.x=unit(c(1), "lines")) +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
      #theme(strip.background = element_blank(), strip.text.x = element_blank())
      ## save plot
      ggsave(file = c(file.path(save_plot_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 4 , width = 16)
    } else {
      D$feature <- factor(D$feature, levels=c("5'UTR", "5' CDS", "3' CDS", "3'UTR"))
      ## plot
      ggplot(D, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity") +
        scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background = element_blank(),
              strip.text.x = element_blank(), text = element_text(size=30),
              legend.justification = c(0.5, 0.9), legend.position = c(0.5, 0.9), legend.direction = "horizontal") +
        guides(fill=guide_legend(title="footprint\nlength")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        facet_grid(~ feature, scales = "free_x", space = "free_x") + scale_x_continuous(breaks= pretty_breaks(2)) +
        theme(panel.spacing.x=unit(c(0.5,1,0.5), "lines")) +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
      #theme(strip.background = element_blank(), strip.text.x = element_blank())
      ## save plot
      ggsave(file = c(file.path(save_plot_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 4 , width = 20)
    }
  }
}




#############################################################################################################
##### stop codon: f.l.d (+ save logo: total and for different stop codons) (for yeast: on CDS only) #####
### geom_bar(stat="identity", position="fill")

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
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/stop_codon"
save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stop"

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
  
  ### logo
  fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
  }
  if (org == "yeast") {
    logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-29):y]}, as.character(txlen$tx_name), txlen$utr5_cds_len))
  } else {
    logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-17):(y+15)]}, as.character(txlen$tx_name), txlen$utr5_cds_len))
  }
  
  logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
  if (org == "zebrafish") {
    logo <- sapply(logo, function(x){paste0(x, collapse="")})
    logo <- logo[nchar(logo) == 33]
  } else {
    logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
  }
  # write to file
  fileConn<-file(c(file.path(save_logo_path, paste0(org, ".txt"))))
  writeLines(logo, fileConn)
  close(fileConn)
  ###
  if (org == "yeast") {
    stop_codons <- sapply(logo, function(x){substr(x, 28, 30)})
  } else {
    stop_codons <- sapply(logo, function(x){substr(x, 16, 18)})
  }
  
  taa <- names(stop_codons[stop_codons == "taa"])
  tga <- names(stop_codons[stop_codons == "tga"])
  tag <- names(stop_codons[stop_codons == "tag"])
 
  
  if (org == "yeast") {
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
    DT_stop <- as.data.frame(gr)
    setDT(DT_stop)
    
    # stop
    DT_stop <- DT_stop[seqnames %in% txlen_stop$tx_name,]
    DT_stop <- DT_stop[grdf_coord_stop,]
    
    # plot per footprint length (from -10 nt if utr5 is available)
    if (org == "yeast") {
      DT_stop$group <- -30:-1
    } else {
      DT_stop$group <- c(-18:15)[-19]
    }
    
    fld_stop <- DT_stop[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stop)[7:(ncol(DT_stop)-1)] ]
    colnames(fld_stop)[2:ncol(fld_stop)] <- sapply(colnames(fld_stop)[2:ncol(fld_stop)], function(x){substr(x,2,3)})
    colourCount <- ncol(fld_stop)-1 # add extra shades to the scale
    fld_stop <- melt(fld_stop, id.var="group")
    
    # add extra shades to the scale
    # colourCount <- ncol(fld_start)-1 # number of levels
    getPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))
    
    ggplot(fld_stop, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_x_continuous(breaks= pretty_breaks(2)) +
      theme(panel.spacing.x=unit(c(1), "lines")) +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ## save plot
    ggsave(file = c(file.path(save_plot_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
    ### make plots for each stop codon separately
    
    ### TAA
    fld_taa <- DT_stop[seqnames %in% taa]
    fld_taa <- fld_taa[, lapply(.SD, sum), by=group, .SDcols=colnames(fld_taa)[7:(ncol(fld_taa)-1)] ]
    colnames(fld_taa)[2:ncol(fld_taa)] <- sapply(colnames(fld_taa)[2:ncol(fld_taa)], function(x){substr(x,2,3)})
    fld_taa <- melt(fld_taa, id.var="group")
    ## plot
    ggplot(fld_taa, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_x_continuous(breaks= pretty_breaks(2)) +
      theme(panel.spacing.x=unit(c(1), "lines")) +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ## save plot
    ggsave(file = c(file.path(save_plot_path, "taa",  paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
    ### TGA
    fld_tga <- DT_stop[seqnames %in% tga]
    fld_tga <- fld_tga[, lapply(.SD, sum), by=group, .SDcols=colnames(fld_tga)[7:(ncol(fld_tga)-1)] ]
    colnames(fld_tga)[2:ncol(fld_tga)] <- sapply(colnames(fld_tga)[2:ncol(fld_tga)], function(x){substr(x,2,3)})
    fld_tga <- melt(fld_tga, id.var="group")
    ## plot
    ggplot(fld_tga, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_x_continuous(breaks= pretty_breaks(2)) +
      theme(panel.spacing.x=unit(c(1), "lines")) +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ## save plot
    ggsave(file = c(file.path(save_plot_path, "tga",  paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
    ### TAG
    fld_tag <- DT_stop[seqnames %in% tag]
    fld_tag <- fld_tag[, lapply(.SD, sum), by=group, .SDcols=colnames(fld_tag)[7:(ncol(fld_tag)-1)] ]
    colnames(fld_tag)[2:ncol(fld_tag)] <- sapply(colnames(fld_tag)[2:ncol(fld_tag)], function(x){substr(x,2,3)})
    fld_tag <- melt(fld_tag, id.var="group")
    ## plot
    ggplot(fld_tag, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(fill=guide_legend(title="footprint\nlength")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_x_continuous(breaks= pretty_breaks(2)) +
      theme(panel.spacing.x=unit(c(1), "lines")) +
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    ## save plot
    ggsave(file = c(file.path(save_plot_path, "tag",  paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)

  }
  
}


##### stall site: f.l.d (+ save logo) #####
### geom_bar(stat="identity", position="fill")
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

paths <- data.frame(consensus = c("/Volumes/USELESS/STALLING/consensus/median/5", "/Volumes/USELESS/STALLING/consensus/median/8",
                                  "/Volumes/USELESS/STALLING/consensus/threshold/5", "/Volumes/USELESS/STALLING/consensus/threshold/8"),
                    plot = c("/Volumes/USELESS/STALLING/PLOTS/stall_site/median/5", "/Volumes/USELESS/STALLING/PLOTS/stall_site/median/8",
                             "/Volumes/USELESS/STALLING/PLOTS/stall_site/threshold/5", "/Volumes/USELESS/STALLING/PLOTS/stall_site/threshold/8"),
                    logo = c("/Volumes/USELESS/STALLING/PLOTS/logo_stall/median/5", "/Volumes/USELESS/STALLING/PLOTS/logo_stall/median/8",
                             "/Volumes/USELESS/STALLING/PLOTS/logo_stall/threshold/5", "/Volumes/USELESS/STALLING/PLOTS/logo_stall/threshold/8"),
                    peaks = c("/Volumes/USELESS/STALLING/peaks/median", "/Volumes/USELESS/STALLING/peaks/median",
                              "/Volumes/USELESS/STALLING/peaks/threshold", "/Volumes/USELESS/STALLING/peaks/threshold"))


for (n in 1:nrow(paths)) {
  consensus_path <- as.character(paths$consensus[n])
  save_plot_path <- as.character(paths$plot[n])
  save_logo_path <- as.character(paths$logo[n])
  peaks_path <- as.character(paths$peaks[n])
  
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
    ## consensus <- consensus[consensus$start < (consensus$utr5_cds_len-11)]
    
    ### logos
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
      if (nrow(consensus) > 0) {
        ## subset on 'seqnames' and 'start' columns from consensus, get 15nt around codon with peak
        peak_idx <- DT[consensus, on=.(seqnames,start), which=TRUE, nomatch=0]
        surr_idx <- lapply(peak_idx, function(x){(x-15):(x+17)})
        ssDT <- DT[unlist(surr_idx),]
        ## get every 33 rows, add up each fragment length
        ssDT$group <- 1:33
        fld <- ssDT[, lapply(.SD, sum, na.rm=TRUE), by=group, .SDcols=colnames(ssDT)[7:(ncol(ssDT)-1)] ]
        ## plot
        colnames(fld)[2:ncol(fld)] <- sapply(colnames(fld)[2:ncol(fld)], function(x){substr(x,2,3)})
        fld$group <- c(-15:17)
        d <- melt(fld, id.var="group")
        ## add extra shades to the scale
        colourCount <- ncol(fld)-1 # number of levels
        getPalette <- colorRampPalette(brewer.pal(11, "RdGy"))
        ## plot
        ggplot(d, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") +
          scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdYlBu"))(colourCount)) + 
          theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
          guides(fill=guide_legend(title="footprint\nlength")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          scale_x_continuous(breaks= pretty_breaks(2)) +
          theme(panel.spacing.x=unit(c(1), "lines")) +
          theme(strip.background = element_blank(), strip.text.x = element_blank())
        ## save plot
        ggsave(file = c(file.path(save_plot_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))))
      }
    }
    
  }
  
}



######################### peaks from every library separately - logos
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))

paths <- data.frame(peaks = c("/Volumes/USELESS/STALLING/peaks/median", "/Volumes/USELESS/STALLING/peaks/threshold"),
                   logo = c("/Volumes/USELESS/STALLING/PLOTS/logo_stall_every/median",
                             "/Volumes/USELESS/STALLING/PLOTS/logo_stall_every/threshold"))


for (n in 1:nrow(paths)) {
  peaks_path <- as.character(paths$peaks[n])
  save_logo_path <- as.character(paths$logo[n])
  
  for (i in 1:length(organisms)) {
    org <- organisms[i]
    fasta_cdna <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
    if (org == "zebrafish") {
      names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
    }
    # load peaks
    libs <- list.files(path = peaks_path, pattern = paste0("^", org))
    
    for (j in 1:length(libs)) {
      load(file = c(file.path(peaks_path, libs[j]))) # peaks
      peaks <- peaks[start > 15]
      
      ### logos - 5
      logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-15):(y+17)]}, as.character(peaks$seqnames), peaks$start))
      logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
      logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
      # write to file
      fileConn<-file(c(file.path(save_logo_path, "5", paste0(substr(libs[j],1,nchar(libs[j])-6), ".txt"))))
      writeLines(logo, fileConn)
      close(fileConn)
      
      ### logos - 8
      peaks <- peaks[zscore > 8]
      logo <- t(mapply(function(x,y){fasta_cdna[[x]][(y-15):(y+17)]}, as.character(peaks$seqnames), peaks$start))
      logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
      logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
      # write to file
      fileConn<-file(c(file.path(save_logo_path, "8", paste0(substr(libs[j],1,nchar(libs[j])-6), ".txt"))))
      writeLines(logo, fileConn)
      close(fileConn)
      
    }
    
  }
  
}


##### SHAPES plots #####

# plot as avg

# example: for mouse_none
fld_stop$avg <- apply(fld_stop[,2:13], 1, function(x){ sum(as.numeric(colnames(fld_stop[,2:13])) * x) / sum(x) })
ggplot(fld_stop, aes(group, avg)) + geom_line() + geom_hline(yintercept=mean(fld_stop$avg), linetype="dashed", color = "red")


#############################################################################################################
##### stop codon: f.l.d as average (for yeast: on CDS only) #####

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
save_plot_path <- "/Volumes/USELESS/STALLING/PLOTS/stop_avg"

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
    DT_stop <- as.data.frame(gr)
    setDT(DT_stop)
    
    # stop
    DT_stop <- DT_stop[seqnames %in% txlen_stop$tx_name,]
    DT_stop <- DT_stop[grdf_coord_stop,]
    
    # plot per footprint length (from -10 nt if utr5 is available)
    if (org == "yeast") {
      DT_stop$group <- -30:-1
    } else {
      DT_stop$group <- c(-18:15)[-19]
    }
    
    fld_stop <- DT_stop[, lapply(.SD, sum), by=group, .SDcols=colnames(DT_stop)[7:(ncol(DT_stop)-1)] ]
    colnames(fld_stop)[2:ncol(fld_stop)] <- sapply(colnames(fld_stop)[2:ncol(fld_stop)], function(x){substr(x,2,3)})
   
    n <- ncol(fld_stop)
    fld_stop$avg <- apply(fld_stop[,2:n], 1, function(x){ sum(as.numeric(colnames(fld_stop[,2:n])) * x) / sum(x) })
    ggplot(fld_stop, aes(group, avg)) + geom_line() + geom_hline(yintercept=mean(fld_stop$avg), linetype="dashed", color = "red")

    ## save plot
    ggsave(file = c(file.path(save_plot_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".png"))), height = 5 , width = 7)
    
  }
  
}



### plot distribution over averages for all nucleotides (codons?) in a transcript (exclude first and last)
### use tails to plot horizontal lines on the stop and stall plots

fld_stop$codon <- rep(c(-5:6)[-6], each=3)
dupa <- fld_stop[, mean(avg), by = c("codon")]
ggplot(dupa, aes(codon, V1)) + geom_line()





