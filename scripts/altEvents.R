# altEvents

library(rtracklayer)
library(data.table)
library(GenomicFeatures)

### altEvents
altEvents <- read.table("/Volumes/USELESS/STALLING/alt_spicing/altEvents.bed")
colnames(altEvents) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand")
setDT(altEvents)

gr <- makeGRangesFromDataFrame(altEvents,
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="chrom",
                               start.field="chromStart",
                               end.field="chromEnd",
                               starts.in.df.are.0based=FALSE)

### CSS
css <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
setDT(css)

### genome
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths <- txLengths[order(txLengths$tx_name),]
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

# CDSs with CSSs only
cds <- cds[names(cds) %in% as.character(css[!is.na(css$h_tx)]$h_tx)]
fo <- findOverlaps(cds, gr)
fo_tx <- sort(names(cds[unique(from(fo))]))

st5 <- c("ENST00000522532", "ENST00000376561", "ENST00000256078", "ENST00000558008", "ENST00000359171",
         "ENST00000368069", "ENST00000336023", "ENST00000223369", "ENST00000272163", "ENST00000341423", "ENST00000560055")
st5_gene <- c("NDUFB9", "CTTN", "KRAS", "PCLAF", "NUP98", "COPA", "TUBA1B", "YKT6", "LBR", "HMGB1", "RPL28")

unique(txLengths[txLengths$tx_name %in% fo_tx,]$gene_id)
ensembl_gene <- read.table("/Volumes/USELESS/STALLING/alt_spicing/ensembl_gene.txt")
colnames(ensembl_gene) <- c("ensembl", "gene")

altCSS <- css[!is.na(css$h_ss)]
altCSS <- altCSS[altCSS$gene %in% ensembl_gene$gene]
altCSS$h_ss <- altCSS$h_ss + 1
altCSS$h_ss_end <- altCSS$h_ss + 2

### find where around these the stall sites are the altEvents
cds_alt <- cds[names(cds) %in% fo_tx]

ss_gr <- makeGRangesFromDataFrame(altCSS,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="h_tx",
                               start.field="h_ss",
                               end.field="h_ss_end",
                               starts.in.df.are.0based=FALSE)
ss_genomic <- mapFromTranscripts(ss_gr, cds)
ss_genomic$seqnames <- as.character(seqnames(ss_gr))

# assign stall site to exon
cds_alt <- unlist(cds_alt)
fo_cds_alt_ss <- findOverlaps(cds_alt, ss_genomic)

cds_alt$ss <- rep(NA, length(cds_alt))
cds_alt$ss[from(fo_cds_alt_ss)] <- start(ss_genomic[to(fo_cds_alt_ss)])

# assign alternative event to exon
fo_cds_alt_altEvent <- findOverlaps(cds_alt, gr)

cds_alt$altEvent <- rep(NA, length(cds_alt))
cds_alt$altEvent[from(fo_cds_alt_altEvent)] <- as.character(gr[to(fo_cds_alt_altEvent)]$name)

###### investigate
cds_alt <- split(cds_alt, names(cds_alt))
cds_alt_gene <- sapply(names(cds_alt), function(x){altCSS[altCSS$h_tx == x]$gene}) # gene name

#dupa <- unlist(cds_alt)
#table(dupa[!is.na(dupa$ss) & !is.na(dupa$altEvent)]$altEvent)

events_next_exon <- c()
same_exon <- 0
for (i in 1:length(cds_alt)) {
  tx <- cds_alt[[i]]
  ss_exon_no <- which(!is.na(cds_alt[[i]]$ss))
  for (j in ss_exon_no) {
    if (j != length(tx)) {
      if (!is.na(tx[j+1]$altEvent)) {
        events_next_exon <- c(events_next_exon, tx[j+1]$altEvent)
        same_exon <- same_exon + sum(j %in% which(!is.na(cds_alt[[i]]$altEvent)))
      }
    }
  }
}

table(events_next_exon)





