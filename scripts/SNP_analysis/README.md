### SNP analysis


Use 0_mm_over_stolsites.R to see analysis
Use SNPs.R to recreate files that were leter processed by the 0_mm_cov.sh commands and 0_mm_over_stolsites.R

human_genomic_granges.Rsave - granges with codons where are stall sites (0-based!!)
download SNPs (common_all_20180418.vcf.gz) from the SNPdb

human_table_bam.Rsave - table with transcripts / and with which BAM (human: h1/h2/h3/h4) libraries contain them, 'tx_ss' is 'tx' in granges data

we analyzed here only h2 library

human.bed

#### processing of clinvar snps/common snps

cut -f1,2 < clinvar_20190513.vcf | sort -k1,1 -k2,2 | uniq | wc -l
cut -d ' ' -f1,2,4,5 < clinvar_20190513.vcf | sort -t ' ' -k1,1 -k2,2 -k3,3 -k4,4 | uniq | wc -l

bcftools filter -R human_renamed.bed


