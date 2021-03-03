# If you want to recreate files in ../../DATA/SNP you can use those commands 
# you will propably have to redownload GRCh38 genome and common_all_20180418.vcf.gz 
# repeat the same for clinvar_20190527.vcf.gz


samtools mpileup -BQ0 -f ../../DATA/SNP/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa ../../DATA/SNP/human_h2.bam  | awk '{print gsub(/\./, "", $5)"\t"$1"\t"$2"\t"$3"\t"$4}' > ../../DATA/SNP/human_h2.counts

samtools mpileup -uf ../../DATA/SNP/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa ../../DATA/SNP/human_h2.bam -d 8000 | bcftools call -mv -Oz -o variants.vcf

bcftools mpileup -d 1000 -f ../../DATA/SNP/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa ../../DATA/SNP/human_h2.bam | ~/Soft/bcftools-1.9/bcftools call -mv -Ov -o human_h2.vcf


bcftools filter -R human_renamed.bed common_all_20180418.vcf.gz -Ov -o ../../DATA/SNP/common_all_filtered.vcf
bcftools filter -R human_renamed.bed human_h2.vcf -Ov -o ../../DATA/SNP/human_h2_fitlered.vcf

