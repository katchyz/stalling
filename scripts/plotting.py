#!/usr/bin/python3

''' Extract fragment lengts and sequences around stall sites '''


def get_fragment_lengths_around_stall_sites(ribo_cov, genes_with_peaks, no_lengths, upstream=15, downstream=15):
	""" Given dictionary with ribosome coverage (per footprint length) and set of genes with peaks,
		returns an array of fragments of different lengths around stall sites """
 
	window_len = upstream + downstream + 3
	meta_stall = np.zeros((no_lengths,window_len)) ### for N different lengths
 
	for tx in genes_with_peaks:
 
		for peak in genes_with_peaks[tx]:
			start = peak-upstream
			end = peak+downstream+3
 
			for l in range(no_lengths):
				if start < 0:
					meta_stall[l][abs(start):] += ribo_cov[tx][1].toarray()[l][0:end]
 
				elif end > len(ribo_cov[tx][1].toarray()[l]):
					meta_stall[l][0:window_len-(end-len(ribo_cov[tx][1].toarray()[l]))] +=
						ribo_cov[tx][1].toarray()[l][start:len(ribo_cov[tx][1].toarray()[l])]
 
				else:
					meta_stall[l] += ribo_cov[tx][1].toarray()[l][start:end]
 
	return meta_stall



def get_sequence_logo_around_stall_sites(fasta, genes_with_peaks, upstream, downstream):
	""" Takes fasta sequences and list of stall sites, and number of nucleotides or amino acids
		upstream and downstream of the stall codon, returns lists of sequences for each stage. """
 
	seq_logo = []
 
	for tx in genes_with_peaks:
 
		if tx in fasta:
			for peak in genes_with_peaks[tx]:
				start = peak - upstream
				end = peak + downstream + 3
 
				if start < 0:
					seq_logo.append('N'*abs(start)+fasta[tx][0:end])
 
				elif end > len(fasta[tx]):
					seq_logo.append(fasta[tx][start:len(fasta[tx])]+'N'*(end-len(fasta[tx])))
 
				else:
					seq_logo.append(fasta[tx][start:end])
 
	return seq_logo



