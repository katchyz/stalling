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



