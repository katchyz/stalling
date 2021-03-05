### Calling peaks and conserved stall sites

#### Workflow for finding conserved stall sites:
    * `workflow_CSS.py` - from BAM, GTF and FASTA to conserved stall sites  

#### Reading and manipulating data:
	* Genomic coordinates:
    	* `read_genome_annotations.py` - reading GTF and BED files
	* Ribosome profiling data:  
    	* `read_data.py` - reading wiggle files and mapping the coverage to genomic annotations  
	* Sequence data:
    	* `read_sequence.py` and `translate.py` - reading FASTA files, translating DNA to protein sequence  
 
#### Calling putative and conserved stall sites:
    * `call_peaks.py` - finding peaks (putative stall sites) on transcripts  
    * `conserved_stall_sites.py` - get conserved stall sites   
    * `old_CSS_biomart.py` (get xml writing and querying) - deprecated, download files directly from BioMart 

#### Saving and plotting data:
    * `save_data.py` - writing peaks to BED file, saving data to and reading from gzip pickle python  
    * `plotting.py` - extract fragment length and sequence data around peaks 
