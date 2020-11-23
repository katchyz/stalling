#!/usr/bin/python3

''' Save data '''

import _pickle as pickle
import gzip


def write_peaks_to_BED_file(genes_with_peaks, bedpath, transcripts, lib='library_name'):
	""" Takes genes_with_peaks, genomic/mRNA coordinates as transcripts and output file path,
		writes peak positions to a BED file. """
 
	bed = open(bedpath, 'w')
	bed.write('track name="Stall sites" visibility=2 itenRgb="On"\n')
 
	for tx in genes_with_peaks:
 
		# stitch the geneToInterval intervals into mRNA with genomic positions
		mRNA = []
		for i in range(len(transcripts[tx]['exon_starts'])):
			mRNA.extend(range(transcripts[tx]['exon_starts'][i], transcripts[tx]['exon_ends'][i]))
 
		CDS_offset = transcripts[tx]['cds_coord'][0] # start codon position
 
		for i in range(len(genes_with_peaks[tx])):
 
			# offset by genes_with_peaks_uniq[stage][tr][i]
			if transcripts[tx]['strand'] == '+':
				peak_pos = CDS_offset + genes_with_peaks[tx][i]
				# get genomic positions
				if len(mRNA) >= peak_pos+3:
					gen_st_pos = mRNA[peak_pos]
					gen_end_pos = mRNA[peak_pos+3]
					bed.write(transcripts[tx]['chromosome']+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+
						lib+'\t0\t'+transcripts[tx]['strand']+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'200,100,0'+'\n')
				else:
					print(f'Peak outside of coordinates of transcript {tx}')
 
			elif transcripts[tx]['strand'] == '-':
				peak_pos = CDS_offset + genes_with_peaks[tx][i] + 3
				if len(mRNA) >= peak_pos:
					gen_st_pos = mRNA[-peak_pos]
					gen_end_pos = mRNA[-peak_pos+3]
					bed.write(transcripts[tx]['chromosome']+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+
						lib+'\t0\t'+transcripts[tx]['strand']+'\t'+str(gen_st_pos)+'\t'+str(gen_end_pos)+'\t'+'0,100,200'+'\n')
				else:
					print(f'Peak outside of coordinates of transcript {tx}')


def save_data(data, filepath=filename.p.gz):
	""" Saves data in python format.
		filepath - filename.p.gz """
 
	if filepath.endswith('p.gz'):
		pickle.dump(data, gzip.open(filepath, 'wb'))
	else:
		print('Provide filename ending with .p.gz')


def load_data(filepath):
	""" Loads data in saved in python gzipped pickle format.
		filepath - filename.p.gz """
 
	if filepath.endswith('p.gz'):
		data = pickle.load(gzip.open(filepath, 'rb'))
	else:
		print('File is not of pickle gzip format. Provide filename ending with .p.gz')



organisms = ["yeast", "fruitfly", "zebrafish", "mouse", "human"]

def export_CSS(conserved_stall_sites, organisms, filename='upsetr.csv'):
	""" Exports conserved stall sites into CSV (for UpSetR plotting) """
 
	f = open(filename, 'w')
	f.write('gene_ss'+','+','.join([str(x) for x in organisms])+'\n')
 
	for gene in conserved_stall_sites:
		for peak_pos in conserved_stall_sites[gene]:
			f.write(gene+'_'+peak_pos+',')
			consorg = conserved_stall_sites[gene][peak_pos].keys()
			#[org in d for org in organisms]
			boo = [int(elem) for elem in [org in consorg for org in organisms]]
			f.write(','.join([str(x) for x in boo]))
			f.write('\n')
 
	f.close()




write_peaks_to_BED_file(gwp, '/Users/kasia/Desktop/test.bed', transcripts, lib='test_bed')



