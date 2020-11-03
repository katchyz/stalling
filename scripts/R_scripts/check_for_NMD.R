### check for peaks on out-of frame stop codons
### (that have a few exons downstream of that stop codon)
library(GenomicFeatures)
library(seqinr)

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

gtf_path <- "/Volumes/USELESS/DATA/genomes/GTF"
org_gtf_fasta <- data.frame(organism = c("yeast", "fruitfly", "zebrafish", "mouse", "human"),
                            gtf = c("Saccharomyces_cerevisiae.R64-1-1.79.gtf", "Drosophila_melanogaster.BDGP6.79.gtf",
                                    "Danio_rerio.GRCz10.81_chr.gtf", "Mus_musculus.GRCm38.79.chr.gtf", "Homo_sapiens.GRCh38.79.chr.gtf"),
                            fasta = c("/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",
                                      "/Volumes/USELESS/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"))


#consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/median/5"
consensus_path <- "/Volumes/USELESS/STALLING/consensus_all/threshold/5"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  load(file = c(file.path(consensus_path, paste0(org, ".Rsave")))) # consensus
  fasta_cds <- read.fasta(as.character(org_gtf_fasta[org_gtf_fasta$organism == org,]$fasta))
  if (org == "zebrafish") {
    names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
  }
  fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]
  fasta_cds <- fasta_cds[names(fasta_cds) %in% unique(as.character(consensus$seqnames))]
  fasta_cds <- lapply(fasta_cds, function(x){paste(x[1:(length(x)-3)], collapse = "")})
  
  # find positions of out-of-frame stop codons (all will be out-of-frame)
  # gregexpr(pattern = 'taa|tga|tag', string)
  stop_codons <- lapply(fasta_cds, function(x){unlist(gregexpr(pattern = 'taa|tga|tag', x))})
  # check if there is a peak within 3nt (before - within 6?)
  
  # ALTERNATIVELY: find stop codons in frames 1 and 2 (not 0)
  # add riboseq for out-of-frame-codons (in frames 1 and 2), check if the z-score is high
  # ...

  
}
