### DISORDER
library(seqinr)
library(GenomicRanges)
library(data.table)
library(GenomicFeatures)

# get DisEMBL results
disembl <- read.fasta(file = "/Volumes/USELESS/STALLING/kmer/css.disembl")

coils <- disembl[grepl("COILS", getAnnot(disembl))]
names(coils) <- sapply(getAnnot(coils), function(x){substring(x, 3, 17)})

rem465 <- disembl[grepl("REM465", getAnnot(disembl))]
names(rem465) <- sapply(getAnnot(rem465), function(x){substring(x, 3, 17)})

hotloops <- disembl[grepl("HOTLOOPS", getAnnot(disembl))]
names(hotloops) <- sapply(getAnnot(hotloops), function(x){substring(x, 3, 17)})

# get intervals
coils_dis <- sapply(getAnnot(coils), function(x){strsplit(substring(x, regexpr("COILS ", x) + 6), split = ", ")})
names(coils_dis) <- names(coils)

rem465_dis <- sapply(getAnnot(rem465), function(x){strsplit(substring(x, regexpr("REM465 ", x) + 7), split = ", ")})
names(rem465_dis) <- names(rem465)

hotloops_dis <- sapply(getAnnot(hotloops), function(x){strsplit(substring(x, regexpr("HOTLOOPS ", x) + 9), split = ", ")})
names(hotloops_dis) <- names(hotloops)

# convert to interval file
coils_dis <- coils_dis[sapply(coils_dis, length) > 0]
coils_start <- sapply(coils_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[1]})})
coils_end <- sapply(coils_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[2]})})
coils_ranges <- IRangesList(mapply(function(x,y){IRanges(as.numeric(x), as.numeric(y))}, coils_start, coils_end))

rem465_dis <- rem465_dis[sapply(rem465_dis, length) > 0]
rem465_start <- sapply(rem465_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[1]})})
rem465_end <- sapply(rem465_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[2]})})
rem465_ranges <- IRangesList(mapply(function(x,y){IRanges(as.numeric(x), as.numeric(y))}, rem465_start, rem465_end))

hotloops_dis <- hotloops_dis[sapply(hotloops_dis, length) > 0]
hotloops_start <- sapply(hotloops_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[1]})})
hotloops_end <- sapply(hotloops_dis, function(x){sapply(strsplit(x, split = "-"), function(y){y[2]})})
hotloops_ranges <- IRangesList(mapply(function(x,y){IRanges(as.numeric(x), as.numeric(y))}, hotloops_start, hotloops_end))

# get CSSs
css <- read.csv(file = "/Volumes/USELESS/STALLING/conservation/ALL_PEAKS/conserved_stall_sites.csv", header = F)
colnames(css) <- c("gene_ss", "y_tx", "y_ss", "f_tx", "f_ss", "z_tx", "z_ss", "m_tx", "m_ss", "h_tx", "h_ss")
css$gene <- sapply(as.character(css$gene_ss), function(x){strsplit(x, split = "_")[[1]][1]})
css <- data.table(tx = as.character(css[!is.na(css$h_ss),]$h_tx), ss = css[!is.na(css$h_ss),]$h_ss + 1)
css$end <- css$ss + 2
css$aa <- floor(css$ss / 3) + 1

# css <- css[css$ss > 45,]

css_ranges <- IRangesList(mapply(function(x,y){IRanges(x,y)}, css$aa, css$aa))
names(css_ranges) <- css$tx


# find overlaps
coils_ranges <- coils_ranges[names(coils_ranges) %in% names(css_ranges)]
coils_css_ranges <- css_ranges[names(css_ranges) %in% names(coils_ranges)]
coils_ranges <- coils_ranges[order(names(coils_ranges))]
coils_css_ranges <- coils_css_ranges[order(names(coils_css_ranges))]
coils_css_ranges <- IRangesList(sapply(unique(names(coils_css_ranges)), function(x){unlist(coils_css_ranges[names(coils_css_ranges) == x])}))
sum(sapply(coils_css_ranges, length)) # 2349 CSS
coils_fo <- mapply(function(x,y){findOverlaps(x,y)}, coils_ranges, coils_css_ranges)
sum(sapply(coils_fo, length)) # 1545 CSSs in coils

rem465_ranges <- rem465_ranges[names(rem465_ranges) %in% names(css_ranges)]
rem465_css_ranges <- css_ranges[names(css_ranges) %in% names(rem465_ranges)]
rem465_ranges <- rem465_ranges[order(names(rem465_ranges))]
rem465_css_ranges <- rem465_css_ranges[order(names(rem465_css_ranges))]
rem465_css_ranges <- IRangesList(sapply(unique(names(rem465_css_ranges)), function(x){unlist(rem465_css_ranges[names(rem465_css_ranges) == x])}))
sum(sapply(rem465_css_ranges, length)) # 2034
rem465_fo <- mapply(function(x,y){findOverlaps(x,y)}, rem465_ranges, rem465_css_ranges)
sum(sapply(rem465_fo, length)) # 331

hotloops_ranges <- hotloops_ranges[names(hotloops_ranges) %in% names(css_ranges)]
hotloops_css_ranges <- css_ranges[names(css_ranges) %in% names(hotloops_ranges)]
hotloops_ranges <- hotloops_ranges[order(names(hotloops_ranges))]
hotloops_css_ranges <- hotloops_css_ranges[order(names(hotloops_css_ranges))]
hotloops_css_ranges <- IRangesList(sapply(unique(names(hotloops_css_ranges)), function(x){unlist(hotloops_css_ranges[names(hotloops_css_ranges) == x])}))
sum(sapply(hotloops_css_ranges, length)) # 2332
hotloops_fo <- mapply(function(x,y){findOverlaps(x,y)}, hotloops_ranges, hotloops_css_ranges)
sum(sapply(hotloops_fo, length)) # 847


tx_coils <- names(coils_fo[sapply(coils_fo, function(x){length(x) > 0})])
tx_rem465 <- names(rem465_fo[sapply(rem465_fo, function(x){length(x) > 0})])
tx_hotloops <- names(hotloops_fo[sapply(hotloops_fo, function(x){length(x) > 0})])
length(unique(c(tx_coils, tx_rem465, tx_hotloops))) # 1339 unique CSSs in disordered regions from coils, rem465 and hotloops

# check control
gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
txdb_h <- makeTxDbFromGFF(c(file.path(gtf_path, "Homo_sapiens.GRCh38.79.chr.gtf")), format="gtf")
txLengths_h <- transcriptLengths(txdb_h, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txLengths_h <- txLengths_h[order(txLengths_h$tx_name),]
css$cds_len <- sapply(css$tx, function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

### CONTROL
utx <- unique(as.character(css$tx))
random_pool <- list()
for (i in 1:length(utx)) {
  txn <- utx[i]
  cds_len_aa <- floor(css[css$tx == txn,]$cds_len[1] / 3)
  around_ss <- c()
  for (ss in css[as.character(css$tx) == txn,]$aa) {
    around_ss <- c(around_ss, c((ss-5):(ss+6)))
  }
  pool <- c(15:(cds_len_aa-15))[!c(15:(cds_len_aa-15)) %in% around_ss]
  random_pool[[txn]] <- pool
}

random_pool <- random_pool[sapply(random_pool, function(x){length(x)}) > 0] # all 1729

css$random1_aa <- sapply(as.character(css$tx), function(x){sample(random_pool[[x]],1)})
css$random2_aa <- sapply(as.character(css$tx), function(x){sample(random_pool[[x]],1)})

random1_ranges <- IRangesList(mapply(function(x,y){IRanges(x,y)}, css$random1_aa, css$random1_aa))
names(random1_ranges) <- css$tx
random2_ranges <- IRangesList(mapply(function(x,y){IRanges(x,y)}, css$random2_aa, css$random2_aa))
names(random2_ranges) <- css$tx

coils_ranges <- coils_ranges[names(coils_ranges) %in% names(random1_ranges)]
coils_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(coils_ranges)]
coils_ranges <- coils_ranges[order(names(coils_ranges))]
coils_random1_ranges <- coils_random1_ranges[order(names(coils_random1_ranges))]
coils_random1_ranges <- IRangesList(sapply(unique(names(coils_random1_ranges)), function(x){unlist(coils_random1_ranges[names(coils_random1_ranges) == x])}))
sum(sapply(coils_random1_ranges, length)) # 2349
coils_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, coils_ranges, coils_random1_ranges)
sum(sapply(coils_fo_random1, length)) # 1339

rem465_ranges <- rem465_ranges[names(rem465_ranges) %in% names(random1_ranges)]
rem465_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(rem465_ranges)]
rem465_ranges <- rem465_ranges[order(names(rem465_ranges))]
rem465_random1_ranges <- rem465_random1_ranges[order(names(rem465_random1_ranges))]
rem465_random1_ranges <- IRangesList(sapply(unique(names(rem465_random1_ranges)), function(x){unlist(rem465_random1_ranges[names(rem465_random1_ranges) == x])}))
sum(sapply(rem465_random1_ranges, length)) # 2034
rem465_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, rem465_ranges, rem465_random1_ranges)
sum(sapply(rem465_fo_random1, length)) # 310

hotloops_ranges <- hotloops_ranges[names(hotloops_ranges) %in% names(random1_ranges)]
hotloops_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(hotloops_ranges)]
hotloops_ranges <- hotloops_ranges[order(names(hotloops_ranges))]
hotloops_random1_ranges <- hotloops_random1_ranges[order(names(hotloops_random1_ranges))]
hotloops_random1_ranges <- IRangesList(sapply(unique(names(hotloops_random1_ranges)), function(x){unlist(hotloops_random1_ranges[names(hotloops_random1_ranges) == x])}))
sum(sapply(hotloops_random1_ranges, length)) # 2332
hotloops_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, hotloops_ranges, hotloops_random1_ranges)
sum(sapply(hotloops_fo_random1, length)) # 624

coils_ranges <- coils_ranges[names(coils_ranges) %in% names(random2_ranges)]
coils_random2_ranges <- random2_ranges[names(random2_ranges) %in% names(coils_ranges)]
coils_ranges <- coils_ranges[order(names(coils_ranges))]
coils_random2_ranges <- coils_random2_ranges[order(names(coils_random2_ranges))]
coils_random2_ranges <- IRangesList(sapply(unique(names(coils_random2_ranges)), function(x){unlist(coils_random2_ranges[names(coils_random2_ranges) == x])}))
sum(sapply(coils_random2_ranges, length)) # 2349
coils_fo_random2 <- mapply(function(x,y){findOverlaps(x,y)}, coils_ranges, coils_random2_ranges)
sum(sapply(coils_fo_random2, length)) # 1254

rem465_ranges <- rem465_ranges[names(rem465_ranges) %in% names(random2_ranges)]
rem465_random2_ranges <- random2_ranges[names(random2_ranges) %in% names(rem465_ranges)]
rem465_ranges <- rem465_ranges[order(names(rem465_ranges))]
rem465_random2_ranges <- rem465_random2_ranges[order(names(rem465_random2_ranges))]
rem465_random2_ranges <- IRangesList(sapply(unique(names(rem465_random2_ranges)), function(x){unlist(rem465_random2_ranges[names(rem465_random2_ranges) == x])}))
sum(sapply(rem465_random2_ranges, length)) # 2034
rem465_fo_random2 <- mapply(function(x,y){findOverlaps(x,y)}, rem465_ranges, rem465_random2_ranges)
sum(sapply(rem465_fo_random2, length)) # 263

hotloops_ranges <- hotloops_ranges[names(hotloops_ranges) %in% names(random2_ranges)]
hotloops_random2_ranges <- random2_ranges[names(random2_ranges) %in% names(hotloops_ranges)]
hotloops_ranges <- hotloops_ranges[order(names(hotloops_ranges))]
hotloops_random2_ranges <- hotloops_random2_ranges[order(names(hotloops_random2_ranges))]
hotloops_random2_ranges <- IRangesList(sapply(unique(names(hotloops_random2_ranges)), function(x){unlist(hotloops_random2_ranges[names(hotloops_random2_ranges) == x])}))
sum(sapply(hotloops_random2_ranges, length)) # 2332
hotloops_fo_random2 <- mapply(function(x,y){findOverlaps(x,y)}, hotloops_ranges, hotloops_random2_ranges)
sum(sapply(hotloops_fo_random2, length)) # 643

###### random

for (i in 1:10) {
  css$random1_aa <- sapply(as.character(css$tx), function(x){sample(random_pool[[x]],1)})
  
  random1_ranges <- IRangesList(mapply(function(x,y){IRanges(x,y)}, css$random1_aa, css$random1_aa))
  names(random1_ranges) <- css$tx
  
  coils_ranges <- coils_ranges[names(coils_ranges) %in% names(random1_ranges)]
  coils_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(coils_ranges)]
  coils_ranges <- coils_ranges[order(names(coils_ranges))]
  coils_random1_ranges <- coils_random1_ranges[order(names(coils_random1_ranges))]
  coils_random1_ranges <- IRangesList(sapply(unique(names(coils_random1_ranges)), function(x){unlist(coils_random1_ranges[names(coils_random1_ranges) == x])}))
  coils_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, coils_ranges, coils_random1_ranges)
  print("coils")
  print(i)
  print(sum(sapply(coils_fo_random1, length)))
  
  rem465_ranges <- rem465_ranges[names(rem465_ranges) %in% names(random1_ranges)]
  rem465_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(rem465_ranges)]
  rem465_ranges <- rem465_ranges[order(names(rem465_ranges))]
  rem465_random1_ranges <- rem465_random1_ranges[order(names(rem465_random1_ranges))]
  rem465_random1_ranges <- IRangesList(sapply(unique(names(rem465_random1_ranges)), function(x){unlist(rem465_random1_ranges[names(rem465_random1_ranges) == x])}))
  rem465_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, rem465_ranges, rem465_random1_ranges)
  print("rem465")
  print(i)
  print(sum(sapply(rem465_fo_random1, length)))
  
  hotloops_ranges <- hotloops_ranges[names(hotloops_ranges) %in% names(random1_ranges)]
  hotloops_random1_ranges <- random1_ranges[names(random1_ranges) %in% names(hotloops_ranges)]
  hotloops_ranges <- hotloops_ranges[order(names(hotloops_ranges))]
  hotloops_random1_ranges <- hotloops_random1_ranges[order(names(hotloops_random1_ranges))]
  hotloops_random1_ranges <- IRangesList(sapply(unique(names(hotloops_random1_ranges)), function(x){unlist(hotloops_random1_ranges[names(hotloops_random1_ranges) == x])}))
  hotloops_fo_random1 <- mapply(function(x,y){findOverlaps(x,y)}, hotloops_ranges, hotloops_random1_ranges)
  print("hotloops")
  print(i)
  print(sum(sapply(hotloops_fo_random1, length)))
  
  tx_coils <- names(coils_fo_random1[sapply(coils_fo_random1, function(x){length(x) > 0})])
  tx_rem465 <- names(rem465_fo_random1[sapply(rem465_fo_random1, function(x){length(x) > 0})])
  tx_hotloops <- names(hotloops_fo_random1[sapply(hotloops_fo_random1, function(x){length(x) > 0})])
  print("total unique tx:")
  print(length(unique(c(tx_coils, tx_rem465, tx_hotloops))))
  
}





# reverse disorder into order
reverse_disorder <- function(dis_tx, cds_len_aa) {
  start <- as.numeric(sapply(strsplit(dis_tx, split = "-"), function(y){y[1]}))
  end <- as.numeric(sapply(strsplit(dis_tx, split = "-"), function(y){y[2]}))
  cds_len <- floor(cds_len_aa/3)-1
  
  order_starts <- c()
  order_ends <- c()
  if (start[1] > 1) {
    order_starts <- c(order_starts, 1)
    order_ends <- c(order_ends, start[1]-1)
  }
  if (length(start) > 1) {
    for (i in 1:(length(end)-1)) {
      order_starts <- c(order_starts, end[i]+1)
      order_ends <- c(order_ends, start[i+1]-1)
    }
  }
  if (end[length(end)] < cds_len) {
    order_starts <- c(order_starts, end[length(end)]+1)
    order_ends <- c(order_ends, cds_len)
  }
  print(start)
  print(end)
  print(cds_len)
  print(order_starts)
  print(order_ends)
  ir <- IRanges(order_starts, order_ends)
  return(ir)
}

coils_cds_len <- sapply(names(coils_dis), function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
rem465_cds_len <- sapply(names(rem465_dis), function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})
hotloops_cds_len <- sapply(names(hotloops_dis), function(x){txLengths_h[txLengths_h$tx_name == x,]$cds_len})

coils_order <- IRangesList(mapply(function(x,y){reverse_disorder(x,y)}, coils_dis, coils_cds_len))
rem465_order <- IRangesList(mapply(function(x,y){reverse_disorder(x,y)}, rem465_dis, rem465_cds_len))
hotloops_order <- IRangesList(mapply(function(x,y){reverse_disorder(x,y)}, hotloops_dis, hotloops_cds_len))

# check with exon boundary
load(file = "/Volumes/USELESS/STALLING/PLOTS/exons/CSSs_position_on_exon.Rsave") # dt_ss
coils_dt <- dt_ss[dt_ss$seqnames.1 %in% names(coils_fo[sapply(coils_fo, length) > 0])]
ggplot(coils_dt, aes(x = position_on_exon)) + geom_histogram() + theme_classic()
rem465_dt <- dt_ss[dt_ss$seqnames.1 %in% names(rem465_fo[sapply(rem465_fo, length) > 0])]
ggplot(rem465_dt, aes(x = position_on_exon)) + geom_histogram() + theme_classic()
hotloops_dt <- dt_ss[dt_ss$seqnames.1 %in% names(hotloops_fo[sapply(hotloops_fo, length) > 0])]
ggplot(hotloops_dt, aes(x = position_on_exon)) + geom_histogram() + theme_classic()
ggsave("/Volumes/USELESS/STALLING/PLOTS/exons/hotoops_position_on_exon.png")
hotloops_dt[hotloops_dt$position_on_exon > 0.85]

ut <- unique(c(tx_coils, tx_rem465, tx_hotloops))
common_dt <- dt_ss[dt_ss$seqnames.1 %in% ut]
ggplot(common_dt, aes(x = position_on_exon)) + geom_histogram() + theme_classic()

tx_coils <- names(coils_css_ranges[sapply(coils_css_ranges, length) > 0])
tx_rem465 <- names(rem465_css_ranges[sapply(rem465_css_ranges, length) > 0])
tx_hotloops <- names(hotloops_css_ranges[sapply(hotloops_css_ranges, length) > 0])
tx <- unique(c(tx_coils, tx_rem465, tx_hotloops))
css[css$tx %in% tx]

# check if stall sites are downstream of ordered regions

#########################################################################
#########################################################################
# (or take predicted values?)

# metaplot around CSSs

# read in CSSs
# get random ranges
## for each css read in disorder table for this tx, extract region around, add to coils, rem465s and hotloops

disorder_path <- "/Volumes/USELESS/STALLING/kmer/disorder"
distx <- as.character(read.table("/Volumes/USELESS/STALLING/kmer/tx.txt")$V1)
css <- css[css$tx %in% distx,]
css <- css[css$aa > 15]
css <- css[css$aa < (floor(css&cds_len/3)-15),]

meta_coils <- data.table(codon = (-15:15))
meta_rem465s <- data.table(codon = (-15:15))
meta_hotloops <- data.table(codon = (-15:15))

random_meta_coils <- data.table(codon = (-15:15))
random_meta_rem465s <- data.table(codon = (-15:15))
random_meta_hotloops <- data.table(codon = (-15:15))

for (i in 1:nrow(css)) {
  tx = css[i,]$tx
  aa = css[i,]$aa
  disembl <- read.table(c(file.path(disorder_path, paste0(tx, ".txt"))))
  coils <- disembl$V1
  rem465s <- disembl$V2
  hotloops <- disembl$V3
  meta_coils <- cbind(meta_coils, coils[(aa-15):(aa+15)])
  meta_rem465s <- cbind(meta_rem465s, rem465s[(aa-15):(aa+15)])
  meta_hotloops <- cbind(meta_hotloops, hotloops[(aa-15):(aa+15)])
  
  rc <- rep(0,31)
  rr <- rep(0,31)
  rh <- rep(0,31)
  # control
  for (j in 1:100) {
    random <- sample(random_pool[[tx]],1)
    rc <- rc + coils[(random-15):(random+15)]
    rr <- rr + rem465s[(random-15):(random+15)]
    rh <- rh + hotloops[(random-15):(random+15)]
  }
  rc <- rc/100
  rr <- rr/100
  rh <- rh/100
  
  random_meta_coils <- cbind(random_meta_coils, rc)
  random_meta_rem465s <- cbind(random_meta_rem465s, rr)
  random_meta_hotloops <- cbind(random_meta_hotloops, rh)
  
  print(i)
}

# complete columns

set(meta_coils, j=which(is.na(colSums(meta_coils))), value=NULL)
set(meta_rem465s, j=which(is.na(colSums(meta_rem465s))), value=NULL)
set(meta_hotloops, j=which(is.na(colSums(meta_hotloops))), value=NULL)
set(random_meta_coils, j=which(is.na(colSums(random_meta_coils))), value=NULL)
set(random_meta_rem465s, j=which(is.na(colSums(random_meta_rem465s))), value=NULL)
set(random_meta_hotloops, j=which(is.na(colSums(random_meta_hotloops))), value=NULL)

# plot
plot_path = "/Volumes/USELESS/STALLING/disorder"

coil <- data.table(codon = c(c(-15:15), c(-15:15)),
                   coil = c(apply(meta_coils[,2:ncol(meta_coils)], 1, mean),
                            apply(random_meta_coils[,2:ncol(random_meta_coils)], 1, mean)),
                   type = c(rep("CSS", 31), rep("control", 31)))

rem465 <- data.table(codon = c(c(-15:15), c(-15:15)),
                     rem465 = c(apply(meta_rem465s[,2:ncol(meta_rem465s)], 1, mean),
                            apply(random_meta_rem465s[,2:ncol(random_meta_rem465s)], 1, mean)),
                   type = c(rep("CSS", 31), rep("control", 31)))

hotloop <- data.table(codon = c(c(-15:15), c(-15:15)),
                      hotloop = c(apply(meta_hotloops[,2:ncol(meta_hotloops)], 1, mean),
                            apply(random_meta_hotloops[,2:ncol(random_meta_hotloops)], 1, mean)),
                   type = c(rep("CSS", 31), rep("control", 31)))

th_coils <- 0.516
th_rem465 <- 0.6
th_hotloops <- 0.1204

ggplot(coil, aes(x=codon, y=coil, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), labels = c("CSS", "control")) +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_hline(yintercept=th_coils, linetype="dashed", color = "red")
ggsave(c(file.path(plot_path, "DisEMBL_coil_th.png")))

ggplot(rem465, aes(x=codon, y=rem465, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), labels = c("CSS", "control")) +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_hline(yintercept=th_rem465, linetype="dashed", color = "red")
ggsave(c(file.path(plot_path, "DisEMBL_rem465_th.png")))

ggplot(hotloop, aes(x=codon, y=hotloop, colour=type)) + geom_line(aes(colour = type, linetype = type)) +
  scale_color_manual(values=c("#615388", "black"), labels = c("CSS", "control")) +
  scale_linetype_manual(values = c("solid", "dashed")) + theme_classic() +
  geom_hline(yintercept=th_hotloops, linetype="dashed", color = "red") + ylim(0.075,0.1205)
ggsave(c(file.path(plot_path, "DisEMBL_hotloop_th.png")))





