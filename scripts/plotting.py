# some plots and reading to BED


def get_fragment_lengths_around_stall_sites(data_dict, genes_with_peaks, no_lengths):
	""" Given dictionary with transcript data and set of genes with peaks,
	returns an array of fragments of different lengths around stall sites """
	meta_stall = {}
	for stage in genes_with_peaks:
		meta_stall[stage] = np.zeros((no_lengths,33)) ### for 7 different lengths
		for tr in genes_with_peaks[stage]:
			for peak in genes_with_peaks[stage][tr]:
				start = peak-15
				end = peak+18
				for l in range(no_lengths):
					if start < 0:
						meta_stall[stage][l][abs(start):] += data_dict[stage][tr][1].toarray()[l][0:end]
					elif end > len(data_dict[stage][tr][1].toarray()[l]):
						meta_stall[stage][l][0:33-(end-len(data_dict[stage][tr][1].toarray()[l]))] += data_dict[stage][tr][1].toarray()[l][start:len(data_dict[stage][tr][1].toarray()[l])]
					else:
						meta_stall[stage][l] += data_dict[stage][tr][1].toarray()[l][start:end]
	return meta_stall


def get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, upstream, downstream):
	""" Gets fasta sequences and list of stall sites, and number of nucleotides or amino acids
	upstream and downstream of the stall codon, returns lists of sequences for each stage.
	Upstream should be a negative integer, downstream a positive """
	seq_logo = {}
	for stage in genes_with_peaks:
		seq_logo[stage] = []
		for tr in genes_with_peaks[stage]:
			if tr in fasta:
				for peak in genes_with_peaks[stage][tr]:
					start = peak + upstream ### vary the length
					end = peak + 3 + downstream ### vary the length
					if start < 0:
						seq_logo[stage].append('N'*abs(start)+fasta[tr][0:end])
					elif end > len(fasta[tr]):
						seq_logo[stage].append(fasta[tr][start:len(fasta[tr])]+'N'*(end-len(fasta[tr])))
					else:
						seq_logo[stage].append(fasta[tr][start:end])
	return seq_logo


# write peaks to BED file
def write_peaks_to_BED_file(genes_with_peaks, filepath, geneToInterval, geneLength, orfs):
	""" Given genes with peaks and output file, writes a BED file for different stages """
	bed = open(filepath, 'w')
	bed.write('track name="Stall sites" visibility=2 itenRgb="On"\n')
	for stage in genes_with_peaks:
		for tr in genes_with_peaks[stage]:
			# stitch the geneToInterval intervals into mRNA with genomic positions
			mRNA = []
			for exon in range(len(geneToInterval[tr])):
				mRNA.extend(range(geneToInterval[tr][exon].start, geneToInterval[tr][exon].end))
			# offset by CDS start position
			if geneToInterval[tr][0].strand == '+':
				CDS_offset = orfs[tr][0][0]
			elif geneToInterval[tr][0].strand == '-':
				CDS_offset = geneLength[tr] - orfs[tr][0][1] + 3
			for i in range(len(genes_with_peaks[stage][tr])):
				# offset by genes_with_peaks_uniq[stage][tr][i]
				if geneToInterval[tr][0].strand == '+':
					peak_pos = CDS_offset + genes_with_peaks[stage][tr][i]
					# get genomic positions
					gen_st_pos = mRNA[peak_pos]
					gen_end_pos = mRNA[peak_pos]+3 # w/o checking for split over exons...
					# print in bed
					bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'200,100,0'+'\n')
				elif geneToInterval[tr][0].strand == '-':
					peak_pos = CDS_offset + genes_with_peaks[stage][tr][i] - 3
					gen_st_pos = mRNA[-peak_pos]-3
					gen_end_pos = mRNA[-peak_pos]
					bed.write(geneToInterval[tr][0].chrom+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+stage+'\t0\t'+geneToInterval[tr][0].strand+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'0,100,200'+'\n')

