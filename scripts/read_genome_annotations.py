#!/usr/bin/python3

''' Reads genome annotations (GTF or BED).
	Returns a dictionary of annotations for each transcript:
	annotations[tx] = { gene_id: gene_id
						chromosome: chrom
						strand: strand
						exon_starts: [exon1_start, exon2_start...]
						exon_ends: [exon1_end, exon2_end...]
						cds_coord: [cds_start, cds_end]
					  }

	Note: BED files do not contain gene_id (for later reduction to one transcript per gene).
'''

##### GTF #####

def parseGTFLine(line):
	if len(line) == 0:
		return {}
	elif not line.find('#') == -1:
		return {}
	result = {}
	input = line.split('\t')
 
	result['chromosome'] = input[0]
	result['input_type'] = input[1]
	result['type'] = input[2]
	result['start'] = int(input[3])
	result['end'] = int(input[4])
	result['direction'] = input[6]
 
	attributes = input[8].split(';')
	for attribute in attributes:
		s = attribute.split()
		if len(s) > 1:
			result[s[0]] = s[1].strip("\"")
 
	return result


def parseGTFFile(gtfFileName):
 
	transcripts = {}
	currentGene = {}
	currentTranscript = {}
 
	gtfFile = open(gtfFileName, 'r')
 
	for i, line in enumerate(gtfFile.readlines()):
 
		entry = parseGTFLine(line)
 
		if not 'transcript_biotype' in entry.keys(): #['transcript_biotype'] == 'protein_coding'
			continue
		if not entry['transcript_biotype'] == 'protein_coding':
			continue
 
		current = {}
 
		geneID = entry['gene_id']
		transcriptID = None
		if not entry['type'] == 'gene':
			transcriptID = entry['transcript_id']
 
		current['gene_id'] = geneID
		current['chromosome'] = entry['chromosome']
		current['forward'] = entry['direction']
		current['RNA_min'] = int(entry['start'])
		current['RNA_max'] = int(entry['end'])
 
		if transcriptID and not transcriptID in transcripts:
			currentTranscript = {}
			currentTranscript['chromosome'] = current['chromosome']
			currentTranscript['strand'] = entry['direction']
			currentTranscript['exon_ends'] = []
			currentTranscript['exon_starts'] = []
			currentTranscript['gene_id'] = geneID
			currentTranscript['CDS_min'] = current['RNA_max']
			currentTranscript['CDS_max'] = current['RNA_min']
			transcripts[transcriptID] = currentTranscript
		elif transcriptID:
			currentTranscript = transcripts[transcriptID]
 
		if entry['type'] == 'exon':
			currentTranscript['exon_ends'].append(current['RNA_max'])
			currentTranscript['exon_starts'].append(current['RNA_min'])
 
		elif entry['type'] == 'start_codon':
			if not 'start_codon' in currentTranscript:
				currentTranscript['start_codon'] = current['RNA_min']
 
		elif entry['type'] == 'stop_codon':
			if not 'stop_codon' in currentTranscript:
				currentTranscript['stop_codon'] = current['RNA_max']
 
		elif entry['type'] == 'CDS':
			currentTranscript['CDS_min'] = min(current['RNA_min'], currentTranscript['CDS_min'])
			currentTranscript['CDS_max'] = max(current['RNA_max'], currentTranscript['CDS_max'])
 
	gtfFile.close()
 
 
	return transcripts


def read_GTF(gtfFileName):
	''' Reads genomic and mRNA coordinates from BED file. '''
 
	transcripts = parseGTFFile(gtfFileName)
 
	for tx in transcripts:
		# sort exon coordinates
		transcripts[tx]['exon_starts'].sort()
		transcripts[tx]['exon_ends'].sort()
 
		# get genomic CDS coordinates
		if 'start_codon' in transcripts[tx].keys() and 'stop_codon' in transcripts[tx].keys():
			transcripts[tx]['cds_gen_coord'] = [transcripts[tx]['start_codon'], transcripts[tx]['stop_codon']]
		elif 'CDS_min' in transcripts[tx].keys() and 'CDS_max' in transcripts[tx].keys():
			transcripts[tx]['cds_gen_coord'] = [transcripts[tx]['CDS_min'], transcripts[tx]['CDS_max']]
		else:
			continue
 
		transcripts[tx]['cds_gen_coord'].sort()
 
		# initialize CDS coordinates
		transcripts[tx]['cds_coord'] = []
 
		# find CDS coordinates
		offset = 0
		for i in range(len(transcripts[tx]['exon_starts'])):
			curr_exon_start = transcripts[tx]['exon_starts'][i]
			curr_exon_end = transcripts[tx]['exon_ends'][i]
			cds_start = transcripts[tx]['cds_gen_coord'][0]
			cds_end = transcripts[tx]['cds_gen_coord'][1]
 
			if curr_exon_start <= cds_start and curr_exon_end >= cds_start:
				transcripts[tx]['cds_coord'].append(int(cds_start - curr_exon_start + offset)) # add CDS start
				offset += curr_exon_end - curr_exon_start
			if curr_exon_start <= cds_end and curr_exon_end >= cds_end:
				transcripts[tx]['cds_coord'].append(int(cds_end - curr_exon_start + offset)) # add CDS end
				continue
			else:
				offset += curr_exon_end - curr_exon_start
 
	return transcripts



##### BED #####


def read_BED(bed_file):
	''' Reads genomic and mRNA coordinates from BED file. '''
 
	bed = open(bed_file, "r")
 
	genomic_coord = {}
	cds_coord = {}
 
	for line in bed:
		gene = line.split("\t")
		chrom = gene[0]
		tx_start = int(gene[1])
		geneStr = gene[3] #"%s:::%s" % (gene[3], gene[1])
		strand = gene[5]
		cds_start = int(gene[6]) - tx_start
		cds_end = int(gene[7]) - tx_start
 
		# Non-coding
		if cds_start == cds_end:
			continue
 
		exonLengths = gene[10].split(",")
		exonStarts = gene[11].split(",")
 
		# BED files occasionally end with a ","
		if gene[10][-1] == ",":
			exonLengths.pop()
		if gene[11][-1] == ",":
			exonStarts.pop()
 
		genomic_coord[geneStr] = {}
		genomic_coord[geneStr]['gene_id'] = ''
		genomic_coord[geneStr]['chromosome'] = chrom
		genomic_coord[geneStr]['strand'] = strand
		genomic_coord[geneStr]['exon_starts'] = []
		genomic_coord[geneStr]['exon_ends'] = []
		genomic_coord[geneStr]['cds_coord'] = []
 
		offset = 0
		for exonStart, exonLen in zip(exonStarts, exonLengths):
			exonStart = int(exonStart)
			exonLen = int(exonLen)
			s = tx_start + exonStart
			e = s + exonLen
 
			# exon intervals
			genomic_coord[geneStr]['exon_starts'].append(s)
			genomic_coord[geneStr]['exon_ends'].append(e)
 
			# mRNA coordinates
			if cds_start >= exonStart and cds_start <= exonStart+exonLen:
				new_cds_start = offset + (cds_start - exonStart)
 
			if cds_end >= exonStart and cds_end <= exonStart+exonLen:
				new_cds_end = offset + (cds_end - exonStart)
 
			offset += exonLen
 
		if new_cds_start > new_cds_end:
			continue # CDS length is negative
 
		genomic_coord[geneStr]['cds_coord'].append(new_cds_start)
		genomic_coord[geneStr]['cds_coord'].append(new_cds_end)
 
	bed.close()
 
	return genomic_coord



#######
# GTF
gtf = '/Volumes/USELESS/DATA/genomes/GTF/Mus_musculus.GRCm38.79.chr.gtf'
annotations = read_GTF(gtf)

# BED
bed = '/Volumes/USELESS/DATA/genomes/BED/Mus_musculus.GRCm38.79.chr.bed'
annotations = read_BED(bed)

annotations['ENSMUST00000152946']
annotations['ENSMUST00000026256']

	