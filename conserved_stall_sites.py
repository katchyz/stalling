### find conserved stall sites in different organisms

import cPickle as pickle
import gzip
import subprocess
from itertools import combinations

# get GWP
# library[stage][transcript] = [peak1, peak2...]

our = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_OUR.p.gz', 'rb'))
giraldez = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_Giraldez.p.gz', 'rb'))
fruitfly = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_fruitfly.p.gz', 'rb'))
mouse = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_mouse.p.gz', 'rb'))
yeast = pickle.load(gzip.open('/Volumes/USELESS/META/genes_with_peaks/GWP_yeast.p.gz', 'rb'))

to_delete = {}
for org in organisms:
	to_delete[org] = {}
	gwps = eval(org)
	for stage in gwps:
		to_delete[org][stage] = []
		for tr in gwps[stage]:
			if not tr in fasta[org]:
				to_delete[org][stage].append(tr)

for org in to_delete:
	for stage in to_delete[org]:
		lib = eval(org)
		for i in to_delete[org][stage]:
			del lib[stage][i]


def write_xml(data, file_out, dataset_name='drerio_gene_ensembl', filter_name='ensembl_transcript_id', attribute_name='ensembl_transcript_id', attribute_name2='uniprot_genename'):
	out_fh = open(file_out, 'w')
	out_fh.write('<?xml version="1.0" encoding="UTF-8"?>\n\
		<!DOCTYPE Query>\n\
		<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >\n\
		\t\t\t\n\
		\t<Dataset name = "'+dataset_name+'" interface = "default" >\n\
		\t\t<Filter name = "'+filter_name+'" value = "'+",".join([(",".join([(key2) for key2 in data[key1]])) for key1 in data])+'"/>\n\
		\t\t<Attribute name = "'+attribute_name+'" />\n\
		\t\t<Attribute name = "'+attribute_name2+'" />\n\
		\t</Dataset>\n\
		</Query>\n')
	out_fh.close()


write_xml(our, '/Volumes/USELESS/META/genes_with_peaks/our_query.xml')
write_xml(giraldez, '/Volumes/USELESS/META/genes_with_peaks/giraldez_query.xml')
write_xml(fruitfly, '/Volumes/USELESS/META/genes_with_peaks/fruitfly_query.xml', dataset_name='dmelanogaster_gene_ensembl', filter_name='flybase_transcript_id', attribute_name='flybase_transcript_id', attribute_name2='flybasename_gene')
write_xml(mouse, '/Volumes/USELESS/META/genes_with_peaks/mouse_query.xml', dataset_name='mmusculus_gene_ensembl')
write_xml(yeast, '/Volumes/USELESS/META/genes_with_peaks/yeast_query.xml', dataset_name='scerevisiae_gene_ensembl', attribute_name2='external_gene_name')

# organisms = ['our', 'giraldez', 'fruitfly', 'mouse', 'yeast']
organisms = ['our', 'giraldez', 'fruitfly', 'mouse', 'yeast']

# add path to GWP ?
def runBiomart(organisms):
	'''Takes a list of organism names, runs BioMart queries (xml files saved before) to get gene names. \
	   Returns a dictionary of organisms `ens_gname[organism][gene_name] = list_of_transcript_ids`. '''
	ens_gname = {}
	for org in organisms:
		ens_gname[org] = {}
		query = '/Volumes/USELESS/META/genes_with_peaks/'+org+'_query.xml'
		biomart = subprocess.Popen(['perl', '/Users/kasia/Documents/PhD/scripts/biomart-perl/scripts/webExample.pl', query],stdout=subprocess.PIPE)
		for line in iter(biomart.stdout.readline,''):
			if len(line.split()) == 2:
				if not line.split()[1].lower() in ens_gname[org]:
					ens_gname[org][line.split()[1].lower()] = [line.split()[0]]
				else:
					ens_gname[org][line.split()[1].lower()].append(line.split()[0])
	return ens_gname

ens_gname = runBiomart(organisms)

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
# fasta[zebrafish] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz')
fasta['our'] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz')
fasta['giraldez'] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz')
fasta['fruitfly'] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.cds.all.fa.gz')
fasta['mouse'] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz')
fasta['yeast'] = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz')
### get our, giraldez...

############
### GET GENES WITH PEAKS!!!, change in line 33
gwpeaks = {}
for org in organisms:
	gwpeaks[org] = {}
	dictionary = eval(org)
	for stage in dictionary:
		for tr in dictionary[stage]:
			if not tr in gwpeaks[org]:
				gwpeaks[org][tr] = []
				gwpeaks[org][tr].extend(dictionary[stage][tr])
			else:
				gwpeaks[org][tr].extend(dictionary[stage][tr])

# get unique peaks and sort...
for org in gwpeaks:
	for tr in gwpeaks[org]:
		gwpeaks[org][tr] = sorted(set(gwpeaks[org][tr]))
# pickle.dump(gwpeaks, gzip.open('/Volumes/USELESS/META/gwpeaks_our_gir_ff_mouse.p.gz', 'wb'))
# pickle.load(gwpeaks, gzip.open('/Volumes/USELESS/META/gwpeaks_our_gir_ff_mouse.p.gz', 'rb'))





ens_gname = runBiomart(organisms)

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
pickle.dump(conSS, gzip.open('/Volumes/USELESS/META/conSS.p.gz', 'wb'))
pickle.dump(lostSS, gzip.open('/Volumes/USELESS/META/lostSS.p.gz', 'wb'))
### get zebrafish in one library perhaps...
### get sequences around stall sites



conserved_stall_sites = {}
for gene in conSS:
	conserved_stall_sites[gene] = {}
	for peak in conSS[gene]:
		conserved_stall_sites[gene][peak] = {}
		for org in conSS[gene][peak]:
			conserved_stall_sites[gene][peak][org[0]] = org[1]

pickle.dump(conserved_stall_sites, gzip.open('/Volumes/USELESS/META/conserved_stall_sites.p.gz', 'wb'))

lost_stall_sites = {}
for gene in lostSS:
	lost_stall_sites[gene] = {}
	for peak in lostSS[gene]:
		lost_stall_sites[gene][peak] = {}
		for org_tr in lostSS[gene][peak]:
			lost_stall_sites[gene][peak][org_tr[0]] = org_tr[1]

pickle.dump(lost_stall_sites, gzip.open('/Volumes/USELESS/META/lost_stall_sites.p.gz', 'wb'))


css = conserved_stall_sites
lss = lost_stall_sites

c = 0
for g in css:
	for p in css[g]:
		if 'our' in css[g][p] and 'giraldez' in css[g][p] and 'fruitfly' in css[g][p] and 'mouse' in css[g][p] and 'yeast' in css[g][p]:
			c += 1




##########
# get sequence
css_logo = open('/Volumes/USELESS/META/css_logo.txt', 'w')
for gene in css:
	for peak in css[gene]:
		for org in css[gene][peak]:
			ss = css[gene][peak][org][0]
			tr = css[gene][peak][org][1]
			logo = fasta[org][tr][ss-15:ss+18]
			if len(logo) == 33:
				css_logo.write(logo+'\n')

css_logo.close()


lss_logo = open('/Volumes/USELESS/META/lss_logo.txt', 'w')
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



