rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

library(GenomicRanges)
library(VariantAnnotation)
library(data.table)
library(rtracklayer)

library(BSgenome.hg38v34)

# all files have to be in common seqlevels style...
# vcf are 1,2,3 not chr1, chr2, chr3

# Count number of unique bases in the x
count_bp <- function(x) {
  x <- granges(x)
  sum(lengths(x))
}

snp_only <- function(x) {
  x[width(x) == 1]
}

# SNPs filtered to overlap human_renamed.bed (stall sites)
all <- snp_only(readVcfAsVRanges("common_all_filtered.vcf"))
clinvar <- snp_only(readVcfAsVRanges("clinvar_filtered.vcf"))

# SNPs filtered to overlap human_renamed_random_SNP.bed - random control stall sites
all_control <- snp_only(readVcfAsVRanges("common_all_filtered_random_SNP.vcf"))
clinvar_control <- snp_only(readVcfAsVRanges("clinvar_filtered_random_SNP.vcf"))

# # this is all cds used to find stall sites
# load("H2cds.Rsave") # 38100
# cds <- unlist(H2cds)
# cds <- reduce(cds)
# seqlevelsStyle(cds) <- seqlevelsStyle(all)
# export.bed(cds, "H2cds_renamed.bed")
cds <- import.bed("H2cds_renamed.bed")

# total number of positions in our search space
count_bp(cds) # 5165354
# how many SNPs fall into this search space?
# bcftools filter -R H2cds_renamed.bed common_all_20180418.vcf.gz -Ov -o common_on_cds.vcf
snps_on_cds <- snp_only(readVcfAsVRanges("common_on_cds.vcf"))
snps_on_cds_clinvar <- snp_only(readVcfAsVRanges("clinvar_on_cds.vcf"))
count_bp(snps_on_cds) # 44147
count_bp(snps_on_cds_clinvar) # 47911
# probability that random base is SNP
p_base_is_snp <- round(count_bp(snps_on_cds) / count_bp(cds), 3) # 0.008546752
p_base_is_snp_clinvar <- round(count_bp(snps_on_cds_clinvar) / count_bp(cds), 3) # 0.009275453
# expected count of SNPs in control
control_ss <- rtracklayer::import.bed("human_renamed_random_SNP.bed")
round(p_base_is_snp * count_bp(control_ss), 0)
round(p_base_is_snp_clinvar * count_bp(control_ss), 0)

# expected count of SNPs in stall sites
conserved_ss <- rtracklayer::import.bed("human_renamed.bed")
p_base_is_snp * count_bp(conserved_ss) # 60.86997
p_base_is_snp_clinvar * count_bp(conserved_ss) # 66.05978

# Count of SNPs on the control
count_bp(all_control)
count_bp(clinvar_control)

# Count of SNPs on the conserved stall sites
count_bp(all)
count_bp(clinvar)

# Count of NON-SYNONYMOUS SNPs
getSeq_ <- function(x) {
  temp <- x
  seqlevelsStyle(temp) <- seqlevelsStyle(BSgenome.hg38v34)
  getSeq(BSgenome.hg38v34, temp)
}

get_count_nsyn_SNP <- function(all_control, control_ss) {
  control_ss_seq <- getSeq_(control_ss)
  control_ss_pep <- as.character(translate(control_ss_seq))
  
  if (any(getSeq_(granges(all_control)) != all_control@ref)) {
    stop("Ref is not equal with reference")
  }
  
  x = 0
  for (i in seq_along(all_control)) {
    temp <- granges(all_control)[i]
    over_ <- control_ss %over% temp
    c_ss <- control_ss[over_][1]
    pos <- abs(start(ranges(c_ss)) - start(temp)) + 1 # 0, 1, 2
    c_ss <- unlist(control_ss_seq[over_])
    c_ss_alt <- c_ss
    if (nchar(all_control@alt[i]) != 1) {
      x = x + 1
    } else {
      c_ss_alt[pos] <- all_control@alt[i]
      if (translate(c_ss) != translate(c_ss_alt)){
        x = x + 1
      } 
    }
  }
  return(x)
}

get_count_nsyn_SNP(all_control, control_ss) # 27
get_count_nsyn_SNP(all, conserved_ss) # 31
get_count_nsyn_SNP(clinvar_control, control_ss) # 35
get_count_nsyn_SNP(clinvar, conserved_ss) # 33


# check whether synonymous frequent to rare changes between control/conserved
# if substitution from frequent to rare it is possible there was little tRNA and that caused stalling
freq <- fread("human_codon_frequency.csv")


get_freq <- function(x) {
  return(freq$Fraction[match(as.character(x), freq$Triplet)])
}
get_frq_syn_SNP <- function(all_control, control_ss) {
  control_ss_seq <- getSeq_(control_ss)
  
  if (any(getSeq_(granges(all_control)) != all_control@ref)) {
    stop("Ref is not equal with reference")
  }
  
  x = c()
  for (i in seq_along(all_control)) {
    temp <- granges(all_control)[i]
    over_ <- suppressWarnings(control_ss %over% temp)
    c_ss <- control_ss[over_][1]
    pos <- abs(start(ranges(c_ss)) - start(temp)) + 1 # 0, 1, 2
    c_ss <- unlist(control_ss_seq[over_])
    c_ss_alt <- c_ss
    if (nchar(all_control@alt[i]) == 1) {
      c_ss_alt[pos] <- all_control@alt[i]
      if (translate(c_ss) == translate(c_ss_alt)){ # only synonymous
        x = rbind(x, c(as.character(c_ss), as.numeric(get_freq(c_ss)), 
                       as.character(c_ss_alt), as.numeric(get_freq(c_ss_alt))))
      }
    }
  }
  return(as.data.table(x))
}

# really random control
# codons that we have into other nonsynonymous codons 
set.seed(42)
get_frq_syn_SNP_random <- function(control_ss) {
  control_ss_seq <- getSeq_(control_ss)
  control_ss_freq <- freq$Fraction[match(as.character(control_ss_seq), freq$Triplet)]
  # changing randomly each sequence to its synonymous variant
  x <- c()
  for (i in seq_along(control_ss_seq)) {
    temp <- unlist(control_ss_seq[i])
    pept <- as.character(translate(temp))
    temp_freq <- freq[freq$`Amino acid` == pept, ]
    temp_freq <- temp_freq[temp_freq$Triplet != as.character(temp), ]
    if (dim(temp_freq)[1] > 0) {
      n <- sample.int(seq_len(dim(temp_freq)[1]))
      x <- c(x, temp_freq$Fraction[n])
    } else {
      x <- c(x, NA)
    }
  }
  ret <- data.table(initial_codon = as.character(control_ss_seq),
                    initial_freq = as.numeric(control_ss_freq),
                    rand_freq = as.numeric(x))
  return(ret[complete.cases(ret), ])
}

all_ctr <- get_frq_syn_SNP(all_control, control_ss)
all_ss <- get_frq_syn_SNP(all, conserved_ss)
round(mean(as.numeric(all_ss$V2)), 3)
round(mean(as.numeric(all_ss$V4)), 3)
expected_change <- get_frq_syn_SNP_random(conserved_ss)
round(mean(expected_change$initial_freq), 3)
round(mean(expected_change$rand_freq), 3)

# small overlap between clinvar and common, yet similar results
sum(clinvar %in% all)
sum(clinvar_control %in% all_control)



clinvar_snps_dt <-
  matrix(c(67, 66, 54, 53),
         nrow = 2,
         dimnames = list(Sample = c("Random", "Conserved"),
                         State = c("Expected", "Found")))
clinvar_snps_dt
fisher.test(clinvar_snps_dt)

common_snps_dt <-
  matrix(c(61, 60, 55, 54),
         nrow = 2,
         dimnames = list(Sample = c("Random", "Conserved"),
                         State = c("Expected", "Found")))
common_snps_dt
fisher.test(common_snps_dt)


# rest is some noise?!
mm <- readVcfAsVRanges("human_h2_fitlered.vcf")
db <- readVcfAsVRanges("clinvar_filtered_random_SNP.vcf")
db2 <- readVcfAsVRanges("clin")

all_control2 <- readVcfAsVRanges("~/Downloads/common_all_filtered_random_SNP.vcf")



# filter while reading
# tabix.file <- TabixFile("common_all_20180418.vcf.gz", yieldSize=10000)
# filterVcf(tabix.file, genome, destination.file,
#           prefilters=prefilters, filters=filters)

stol <- rtracklayer::import("human_random_SNP.bed") #"human.bed"
stol <- stol[width(stol) == 3]
rtracklayer::export.bed(stol, "human2.bed")
stol <- rtracklayer::import.bed("human2.bed")

seqlevelsStyle(stol) <- seqlevelsStyle(db)
rtracklayer::export.bed(stol, "human_renamed_random_SNP.bed")
stol <- rtracklayer::import.bed("human_renamed_random_SNP.bed")

mm_over <- mm[mm %over% stol]
mm_over@altDepth/(mm_over@refDepth + mm_over@altDepth)

#load("human_genomic_granges.Rsave")
load("human_genomic_granges.Rsave")
