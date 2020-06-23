### control logos for each library (ligation biases?)

# load peaks
# get positions ~50 nt up/downstream from stall site (check if still within CDS)
# save logos

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))
peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
save_logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/random"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  ### longest tx
  txdb <- makeTxDbFromGFF(c(file.path(gtf_path, as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$gtf))), format="gtf")
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  # strand
  strand <- transcripts(txdb)
  strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
  txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  txLengths <- txLengths[order(rownames(txLengths)),]
  ## longest tx
  txlen <- arrange(txLengths, gene_id, desc(tx_len))
  txlen <- txlen[!duplicated(txlen$gene_id),]
  rownames(txlen) <- txlen$tx_name
  ### fasta
  fasta_cds <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
  }
  # libs
  libs <- list.files(path = peaks_path, pattern = paste0("^", org))
  for (j in 1:length(libs)) {
    load(file = c(file.path(peaks_path, libs[j]))) # peaks
    peaks <- peaks[peaks$seqnames %in% txlen$tx_name,]
    peaks <- peaks[,1:2]
    # add utr5 and cds lengths
    peaks$utr5_len <- sapply(peaks$seqnames, function(x){txLengths[txLengths$tx_name == x,]$utr5_len})
    peaks$cds_len <- sapply(peaks$seqnames, function(x){txLengths[txLengths$tx_name == x,]$cds_len})
    peaks$ss_cds <- peaks$start - peaks$utr5_len
    peaks$random_cds <- peaks$ss_cds + 50
    peaks <- peaks[peaks$random_cds < (peaks$utr5_len + peaks$cds_len),]
    ### extract logo
    logo <- t(mapply(function(x,y){fasta_cds[[x]][(y-15):(y+17)]}, as.character(peaks$seqnames), peaks$random_cds))
    logo <- logo[apply(logo, 1, function(x){!any(is.na(x))}),]
    logo <- apply(logo, 1, function(x){paste0(x, collapse="")})
    # write to file
    fileConn<-file(c(file.path(save_logo_path, paste0(substr(libs[j],1,nchar(libs[j])-6), ".txt"))))
    writeLines(logo, fileConn)
    close(fileConn)
  }
}



