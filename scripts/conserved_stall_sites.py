#!/usr/bin/python3

### find conserved stall sites in different organisms

import os.path
import gzip

import cPickle as pickle
import subprocess
from itertools import combinations

# get GWP
# library[stage][transcript] = [peak1, peak2...]
#path = '/Volumes/USELESS/STALLING/conservation'
path = '/export/valenfs/projects/kasia/STALLING'
organisms = ['yeast', 'fruitfly', 'zebrafish', 'mouse', 'human']
# os.path.join('/my/root/directory', 'in', 'here')


gwpeaks = {}
for org in organisms:
	f = open(os.path.join(path, 'consensus_df/all_threshold', org+'.csv'), 'r')
	gwpeaks[org] = {}
	for line in f:
		tx = eval(line.split(',')[1])
		start = eval(line.split(',')[5])
		if not tx == 'seqnames':
			if not tx in gwpeaks[org].keys():
				gwpeaks[org][tx] = [start]
			else:
				gwpeaks[org][tx].append(start)




def readBiomart(organisms):
	'''Takes a list of organism names, reads BioMart CSV files to get gene names. \
	   Returns a dictionary of organisms `ens_gname[organism][gene_name] = list_of_transcript_ids`. '''
	ens_gname = {}
	for org in organisms:
		ens_gname[org] = {}
		export = open('/export/valenfs/projects/kasia/STALLING/biomart_export/'+org+'.txt', 'r')
		for line in export:
			line = line.strip()
			tx = line.split(',')[0]
			gene = line.split(',')[1].upper()
			if not tx == 'Transcript stable ID':
				if tx in gwpeaks[org].keys():
					if not gene in ens_gname[org]:
						ens_gname[org][gene] = [tx]
					else:
						ens_gname[org][gene].append(tx)
	return ens_gname


ens_gname = readBiomart(organisms[1:])
ens_gname['yeast'] = {}
export = open('/export/valenfs/projects/kasia/STALLING/biomart_export/yeast.txt', 'r')
for line in export:
	line = line.strip()
	tx = line.split(',')[0]
	gene = line.split(',')[1]
	if not gene == 'n/a':
		if not gene == '':
			if tx in gwpeaks['yeast'].keys():
				if not gene in ens_gname['yeast']:
					ens_gname['yeast'][gene] = [tx]
				else:
					ens_gname['yeast'][gene].append(tx)



def parse_clustalo(clustalo_subprocess):
	""" Writes output of clustalo into dictionary """
	alignment = {}
	header_cds = ''
	seq_cds = ''
	for line in iter(clustalo_subprocess.stdout.readline, ''):
		line = line.strip()
		if line.startswith('>'):
			if header_cds and seq_cds:
				alignment[header_cds] = seq_cds
			header_cds = line.split()[0][1:]
			seq_cds = ''
		else:
			seq_cds += line
		# for the last sequence:
		alignment[header_cds] = seq_cds
	return alignment


def find_absolute_peak_position(sequence, peak):
	'''Finds position of a peak in an alignment'''
	peak = peak + 1 # to account for 0-based python indexing
	gaps = 0
	nt = 0
	for i in sequence:
		if i == '-':
			gaps += 1
		else:
			nt += 1
		if nt == peak:
			break
	abs_peak = peak + gaps - 1
	return abs_peak


def get_FASTA_sequence(filepath):
	""" Reads gzipped fasta file """
	f = gzip.open(filepath, 'rb')
	fasta = {}
	header_cds = ''
	seq_cds = ''
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if header_cds and seq_cds:
				fasta[header_cds] = seq_cds
			header_cds = line.split()[0][1:]
			seq_cds = ''
		else:
			seq_cds += line
		# for the last sequence:
		fasta[header_cds] = seq_cds
	return fasta


fasta = {}
fasta['yeast'] = get_FASTA_sequence('/Home/ii/katchyz/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz')
fasta['fruitfly'] = get_FASTA_sequence('/Home/ii/katchyz/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz')
fasta['zebrafish'] = get_FASTA_sequence('/Home/ii/katchyz/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz')
# names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
for key in fasta['zebrafish'].keys():
	fasta['zebrafish'][key[0:18]] = fasta['zebrafish'][key]
	del fasta['zebrafish'][key]

fasta['mouse'] = get_FASTA_sequence('/Home/ii/katchyz/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')
fasta['human'] = get_FASTA_sequence('/Home/ii/katchyz/DATA/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz')





conSS = {} ##### !!!!!!!!!!
lostSS = {}

for n in range(len(organisms), 1, -1):
	for curr_org_list in combinations(organisms, n):
		print 'gene present in: ', curr_org_list
		### GET GENES IN COMMON
		gene_tr_list = [ens_gname[org] for org in curr_org_list] # is it in organisms or curr_org_list?
		genes_in_common = set.intersection(*map(set, gene_tr_list)) # set
		for gene in genes_in_common:
			if not gene in conSS: ##### !!!!!!!!!!
				conSS[gene] = {}
				lostSS[gene] = {}
			### get longest transcript
			### write to homologs.fasta
			homologs = open('homologs.fasta', 'w')
			for org in curr_org_list:
				fasta_list = []
				for i in range(len(ens_gname[org][gene])):
					if ens_gname[org][gene][i] in fasta[org]:
						fasta_list.append(fasta[org][ens_gname[org][gene][i]])
					else:
						fasta_list.append('X') ### if fasta_list == ['X'], do not use this!!!
				max_index = fasta_list.index(max(fasta_list, key=len))
				homologs.write('>'+org+'_'+ens_gname[org][gene][max_index]+'\n'+fasta[org][ens_gname[org][gene][max_index]]+'\n')
			homologs.close()
			### run ClustalO
			clustalo = subprocess.Popen(['clustalo', '-i', 'homologs.fasta'], stdout=subprocess.PIPE)
			alignment = parse_clustalo(clustalo)
			### get absolute and relative peaks
			abs_peaks = {}
			rel_peaks = {}
			for key in alignment:
				org, tr = key.split('_')[0], key.split('_')[1]
				abs_peaks[org] = []
				rel_peaks[org] = []
				for peak in gwpeaks[org][tr]:
					abs_peaks[org].append(find_absolute_peak_position(alignment[key], peak))
					rel_peaks[org].append((peak, tr))
			### GET STALL SITES IN COMMON
			for n2 in range(len(curr_org_list), 1, -1):
				for curr_org_list2 in combinations(curr_org_list, n2):
					print 'stall site conserved in: ', curr_org_list2
					############# GET CODE FOR PEAKS FROM BELOW HERE !!!!
					abs_peaks_list = [(abs_peaks[org], org) for org in curr_org_list2]
					for peak in abs_peaks_list[0][0]:
						indices = {}
						indices[abs_peaks_list[0][1]] = abs_peaks_list[0][0].index(peak)
						conSS[gene][str(peak)] = []
						conSS[gene][str(peak)].append([abs_peaks_list[0][1], rel_peaks[abs_peaks_list[0][1]][indices[abs_peaks_list[0][1]]]]) ##### !!!!!!!!!!
						lostSS[gene][str(peak)] = []
						for i in range(1, len(abs_peaks_list)):
							ps = abs_peaks_list[i][0]
							if ps:
								if peak-3 <= min(ps, key=lambda x:abs(x-peak)) <= peak+3:
									indices[abs_peaks_list[i][1]] = abs_peaks[abs_peaks_list[i][1]].index(min(ps, key=lambda x:abs(x-peak)))
									conSS[gene][str(peak)].append([abs_peaks_list[i][1], rel_peaks[abs_peaks_list[i][1]][indices[abs_peaks_list[i][1]]]]) ##### !!!!!!!!!!
						if len(indices) == len(abs_peaks_list):
							print 'have a match:', peak # + DEL PEAK
							print [(org, rel_peaks[org][indices[org]][0]) for org in indices]
							for org in abs_peaks:
								if org in indices:
									del abs_peaks[org][indices[org]]
									del rel_peaks[org][indices[org]]
						else:
							del conSS[gene][str(peak)]
					### push the org with missing stall site to separate dictionary
						for o in curr_org_list:
							if not o in curr_org_list2 and len(indices) == len(abs_peaks_list): #### AND IF REMAINING ORGANISMS HAVE THE STALL SITE CONSERVED!
								if rel_peaks[o]:
									k = str(o)+'_'+str(rel_peaks[o][0][1])
									seq = alignment[k][0:peak]
									rel_p = len(seq) - seq.count('-')
									lostSS[gene][str(peak)].append([o, (rel_p, rel_peaks[o][0][1])])
						if not lostSS[gene][str(peak)]:
							del lostSS[gene][str(peak)]
			for org in ens_gname:
				if gene in ens_gname[org]:
					del ens_gname[org][gene] ### DELETE
			if not conSS[gene]:
				del conSS[gene]
			if not lostSS[gene]:
				del lostSS[gene]


###
pickle.dump(conSS, gzip.open('/export/valenfs/projects/kasia/STALLING/all_threshold/conSS.p.gz', 'wb'))
pickle.dump(lostSS, gzip.open('/export/valenfs/projects/kasia/STALLING/all_threshold/lostSS.p.gz', 'wb'))


### get sequences around stall sites

conserved_stall_sites = {}
for gene in conSS:
	conserved_stall_sites[gene] = {}
	for peak in conSS[gene]:
		conserved_stall_sites[gene][peak] = {}
		for org in conSS[gene][peak]:
			conserved_stall_sites[gene][peak][org[0]] = org[1]

pickle.dump(conserved_stall_sites, gzip.open('/export/valenfs/projects/kasia/STALLING/all_threshold/conserved_stall_sites.p.gz', 'wb'))

lost_stall_sites = {}
for gene in lostSS:
	lost_stall_sites[gene] = {}
	for peak in lostSS[gene]:
		lost_stall_sites[gene][peak] = {}
		for org_tr in lostSS[gene][peak]:
			lost_stall_sites[gene][peak][org_tr[0]] = org_tr[1]

pickle.dump(lost_stall_sites, gzip.open('/export/valenfs/projects/kasia/STALLING/all_threshold/lost_stall_sites.p.gz', 'wb'))


css = conserved_stall_sites
lss = lost_stall_sites

c = 0
for g in css:
	for p in css[g]:
		if 'yeast' in css[g][p] and 'fruitfly' in css[g][p] and 'zebrafish' in css[g][p] and 'mouse' in css[g][p] and 'human' in css[g][p]:
			c += 1




##########
# get sequence
css_logo = open('/export/valenfs/projects/kasia/STALLING/all_threshold/css_logo.txt', 'w')
for gene in css:
	for peak in css[gene]:
		for org in css[gene][peak]:
			ss = css[gene][peak][org][0]
			tr = css[gene][peak][org][1]
			logo = fasta[org][tr][ss-15:ss+18]
			if len(logo) == 33:
				css_logo.write(logo+'\n')

css_logo.close()


lss_logo = open('/export/valenfs/projects/kasia/STALLING/all_threshold/lss_logo.txt', 'w')
for gene in lss:
	for peak in lss[gene]:
		for org in lss[gene][peak]:
			ss = lss[gene][peak][org][0]
			tr = lss[gene][peak][org][1]
			logo = fasta[org][tr][ss-15:ss+18]
			if len(logo) == 33:
				lss_logo.write(logo+'\n')

lss_logo.close()
## do it on absolute alignment - see if there's any deletions etc.



