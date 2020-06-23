### align colourspace reads (Solid)

bowtie-build -C #...

#tophat -G /Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf --bowtie1 --transcriptome-index=tophat_tx_colour_index bowtie_dna_colour_index


tophat -C -G /Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf --bowtie1 --transcriptome-index=tophat_tx_C_index bowtie_dna_colour_index


### 1cell
bowtie_index='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/bowtie_dna_colour_index'
tophat_trans='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/tophat_tx_C_index/Danio_rerio.GRCz10.84'
raw="/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish"
processed="/export/valenfs/data/processed_data/RNA-seq/aanes_2011_zebrafish"
out=${processed}/aligned_with_colorspace_index
csfa=${raw}/SRR062658_F3.csfasta
qual=${raw}/SRR062658_F3_QV.qual

tophat -C --quals --no-coverage-search --bowtie1 --transcriptome-index=$tophat_trans -o $out --bowtie1 $bowtie_index $csfa $qual




### 3.5h
bowtie_index='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/bowtie_dna_colour_index'
tophat_trans='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/tophat_tx_C_index/Danio_rerio.GRCz10.84'
raw="/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish"
processed="/export/valenfs/data/processed_data/RNA-seq/aanes_2011_zebrafish"
out=${processed}/aligned_with_colorspace_index_3_5h
csfa=${raw}/SRR062661_F3.csfasta
qual=${raw}/SRR062661_F3_QV.qual

tophat -C --quals --no-coverage-search --bowtie1 --transcriptome-index=$tophat_trans -o $out --bowtie1 $bowtie_index $csfa $qual


########### vesterlund, 1 cell

bowtie_index='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/bowtie_dna_colour_index'
tophat_trans='/export/valenfs/data/raw_data/RNA-Seq/aanes_2011_zebrafish/tophat_tx_C_index/Danio_rerio.GRCz10.84'

raw="/export/valenfs/data/raw_data/RNA-Seq/vesterlund_2011_zebrafish"
processed="/export/valenfs/data/processed_data/RNA-seq/vesterlund_2011_zebrafish"
out=${processed}/aligned_with_colorspace_index
csfa=${raw}/KI-BN-JKE-DRERIO-RNASEQ-20090709-1cell_F3.csfasta
qual=${raw}/KI-BN-JKE-DRERIO-RNASEQ-20090709-1cell_F3_QV.qual

tophat -C --quals --no-coverage-search --bowtie1 --transcriptome-index=$tophat_trans -o $out --bowtie1 $bowtie_index $csfa $qual


