#!/usr/bin/python3

''' run_biomart (deprecated, download files from BioMart directly)
	uses biomart-perl script webExample.pl '''


import os.path
import subprocess

path = '../DATA'
organisms = ['yeast', 'fruitfly', 'zebrafish', 'mouse', 'human']



def write_xml(data, file_out, dataset_name='drerio_gene_ensembl', filter_name='ensembl_transcript_id', attribute_name='ensembl_transcript_id', attribute_name2='uniprot_genename'):
	out_fh = open(file_out, 'w')
	out_fh.write('<?xml version="1.0" encoding="UTF-8"?>\n\
		<!DOCTYPE Query>\n\
		<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >\n\
		\t\t\t\n\
		\t<Dataset name = "'+dataset_name+'" interface = "default" >\n\
		\t\t<Filter name = "'+filter_name+'" value = "'+",".join(data.keys())+'"/>\n\
		\t\t<Attribute name = "'+attribute_name+'" />\n\
		\t\t<Attribute name = "'+attribute_name2+'" />\n\
		\t</Dataset>\n\
		</Query>\n')
	out_fh.close()



write_xml(gwpeaks['yeast'], os.path.join(path, 'biomart_queries', 'yeast.xml'), dataset_name='scerevisiae_gene_ensembl', attribute_name2='external_gene_name')
write_xml(gwpeaks['fruitfly'], os.path.join(path, 'biomart_queries', 'fruitfly.xml'), dataset_name='dmelanogaster_gene_ensembl', filter_name='flybase_transcript_id', attribute_name='flybase_transcript_id', attribute_name2='flybasename_gene')
write_xml(gwpeaks['zebrafish'], os.path.join(path, 'biomart_queries', 'zebrafish.xml'))
write_xml(gwpeaks['mouse'], os.path.join(path, 'biomart_queries', 'mouse.xml'), dataset_name='mmusculus_gene_ensembl')
write_xml(gwpeaks['human'], os.path.join(path, 'biomart_queries', 'human.xml'), dataset_name='hsapiens_gene_ensembl')


# add path to GWP ?
def runBiomart(organisms):
	'''Takes a list of organism names, runs BioMart queries (xml files saved before) to get gene names. \
	   Returns a dictionary of organisms `ens_gname[organism][gene_name] = list_of_transcript_ids`. '''
	ens_gname = {}
	for org in organisms:
		ens_gname[org] = {}
		query = os.path.join(path, 'biomart_queries', org, '.xml')
		biomart = subprocess.Popen(['perl', os.path.join(path, 'biomart-perl/scripts/webExample.pl'), query],stdout=subprocess.PIPE)
		for line in iter(biomart.stdout.readline,''):
			if len(line.split()) == 2:
				if not line.split()[1].lower() in ens_gname[org]:
					ens_gname[org][line.split()[1].lower()] = [line.split()[0]]
				else:
					ens_gname[org][line.split()[1].lower()].append(line.split()[0])
	return ens_gname


ens_gname = runBiomart(organisms)



