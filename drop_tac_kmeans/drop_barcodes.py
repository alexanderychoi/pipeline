# -*-coding:Utf-8 -*
## TITAN WILL LIVE FOR EVER IN OUR HEARTS
#######################################
## drop_barcodes is designed to align 
## the plasmid barcodes only, while 
## using the barcode clusters assigned
## with the rest of the transcriptome
## data 
#######################################
## Developed by:
## Sisi Chen 
## sisi.chen1@gmail.com
#######################################

from drop_titan_params_sisi import *
from drop_barcodes_params import * # needs to be in this order to override some parameters from above file
from itertools import izip
import os
import os.path
import time
import gzip
import operator
import glob
import pdb
import Levenshtein
import pickle
from collections import defaultdict


def str_dist(s1, s2):
	# dynamic programming algorithm to compute the 
	# edit distance (levenshtein distance) between two strings
	return Levenshtein.distance(s1,s2)

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

class fastq_read:
	line1 = ''
	line2 = ''
	line3 = ''
	line4 = ''
	umi = ''
	pre_barcode = ''
	post_barcode = ''
	tac_length = len(str_search)
	def __init__(self, line1, line2, line3, line4, f1_line2):
		self.line1 = line1
		self.line2 = line2
		self.line3 = line3
		self.line4 = line4
		self.umi = f1_line2[tac_length+barcode_length:tac_length+barcode_length+umi_length]
		self.pre_barcode = f1_line2[tac_length:tac_length+barcode_length]
	def clear_read(self):
		self.line1 = ''
		self.line2 = ''
		self.line3 = ''
		self.line4 = ''
		self.umi = ''
		self.pre_barcode = ''
		self.post_barcode = ''		
	def remove_polyA(self):
		polyA_seq = 'A'*A_num 
		if polyA_seq in self.line2:
			A_start_ind = str.find(self.line2, polyA_seq)
			self.line2 = self.line2[0:A_start_ind] + '\n'
			self.line4 = self.line4[0:A_start_ind] + '\n'
		self.check_length()			
	def remove_tso(self):
		tso_match = longest_common_substring(tso, self.line2)
		tso_match_start = tso_match[0]
		tso_match_len = tso_match[1]
		if tso_match_len>tso_word_len:
			self.clear_read()
	def check_length(self):
		if len(self.line2)<min_read_len:
			self.clear_read()
	def write_read(self, outfile):
		if self.line1:
			outfile.write(self.line1 + self.line2 + self.line3 + self.line4)

def preprocess_barcodes(f1, f2):
	print "first pass filtering of fastq files:"
	file1_fastq = gzip.open(f1,'rb')
	file2_fastq = gzip.open(f2,'rb')
	print f1.split('/')[-1]
	print f2.split('/')[-1]
	file_processed = open(f1_base_name+'_processed.fastq', 'w+', 1000)
	file_umi = open(f1_base_name+'_umi.txt', 'w+',1)
	file_pre_barcode = open(f1_base_name+'_pre_barcode.txt', 'w+', 1000)
	# #Stats about the trimming process
	total_reads = 0
	f1_line1 = 1
	while f1_line1:
		f1_line1 = file1_fastq.readline()
		f1_line2 = file1_fastq.readline()
		f1_line3 = file1_fastq.readline()
		f1_line4 = file1_fastq.readline()
		f2_line1 = file2_fastq.readline()
		f2_line2 = file2_fastq.readline()
		f2_line3 = file2_fastq.readline()
		f2_line4 = file2_fastq.readline()
		total_reads+=1
		const_region = f2_line2[11:19]
		c_dist = str_dist(const_region, p_str)
		if (c_dist<1):
			curr_read = fastq_read(f2_line1, f2_line2, f2_line3, f2_line4, f1_line2)
			curr_read.remove_polyA()
			if curr_read.line1:
				file_umi.write(curr_read.umi+'\n')
				umi_list.append(curr_read.umi)
				pre_barcode_list.append(curr_read.pre_barcode)			
				file_pre_barcode.write(curr_read.pre_barcode+'\n')				
				curr_read.write_read(file_processed)
	file_pre_barcode.close()
	file_umi.close()		
	file1_fastq.close()
	file2_fastq.close()
	file_processed.close()
	return [umi_list, pre_barcode_list, total_reads]

def get_gene_names_from_fasta(fasta_path):
	# Get the list of all genes from the genome fasta file
	gene_counter = 1 #start coutner at 1 to offset for header
	dict_gene_names = defaultdict(str)
	dict_gene_counter = defaultdict(int)
	fasta_file = open(fasta_path,'r')
	while True:
		line = fasta_file.readline()
		if not line:
			break
		else:
			if line[:1] == '>':
				gene_symbol=line.split(' ')[1]
				gene_symbol=gene_symbol[5:]
				gene_symbol=gene_symbol.strip().upper()
				if gene_symbol not in dict_gene_counter:
					dict_gene_counter[gene_symbol] = gene_counter
					gene_counter+=1
				gene_nm=line.split(' ')[0]
				gene_nm=gene_nm[1:]
				gene_nm=gene_nm.strip().upper()
				dict_gene_names[gene_nm] = gene_symbol
	fasta_file.close()
	return [dict_gene_names, dict_gene_counter]

def align_reads(fastq_file):
	start_time = time.time()
	fastq_file_name = fastq_file.split('/')[-1]
	alignment_file_name = dir_path_alignment + fastq_file.split('/')[-1].split(os.extsep)[0] + '.sam'
	if os.path.isfile(alignment_file_name+'.gz'):
		print "Skipping alignment: " + fastq_file.split(os.extsep)[0] + '.sam.gz already exists...'
	else:
		print "Starting alignment with bowtie2 for", fastq_file_name
		os.system("nice " + bowtie2_dir + " bowtie2 " + bowtie_opt +\
			" -x " + reference_genome + " -U " + fastq_file +\
			" -S " + alignment_file_name)
		print 
		os.system("gzip " + alignment_file_name)
		print fastq_file_name + " aligned and file compressed"
		print "..................................................................."
	total_time = time.time() - start_time
	print "Reads alignment time:"
	print int(total_time/60),"min",int(total_time%60),"sec"
	return alignment_file_name

def txt_file_to_list(filename):
	in_file = open(filename, 'r')#,1000)
	row_list = in_file.read().splitlines()
	in_file.close()
	return row_list

def read_sam_umi(sam_file_name, f1_base_name):

	############################### Data gathering ###############################
	print "**********************************"
	print "**       Retrieving data        **"
	print "**********************************"
	print "\n"

	if os.path.isfile(dir_path_barcodes+'matrix.txt'):
		print "matrix.txt already exists.......................................... 100 %"
	else:

		sam_file = gzip.open(sam_file_name + '.gz', 'rb')
		corrected_bcs = txt_file_to_list(f1_base_name+'_barcode.txt')	
		umi_list = txt_file_to_list(f1_base_name + '_umi.txt')
		list_f = os.listdir(dir_path_fastqs)
		
		barcode_counter = 1
		dict_genes_barcode = defaultdict(lambda : defaultdict(int))	
		dict_barcode_counter = defaultdict(int)
		dict_barcode_occurences = defaultdict(int)

		print "Storing data in dictionaries..."

		# collect alignment statistics
		sam_gene=0
		sam_star=0
		bowtie_al=0
		dict_quality = defaultdict(lambda : defaultdict(int))
		dict_quality_scores=defaultdict(int)
		sam_line_ind = 0
		dis_redund = 0
		read_has_no_bc = 0

		# map gene to tuple containing all umi and barcode pairs
		gene_to_umi_bc_dict = defaultdict(lambda: dict())

		while True: #sam_line_ind<1000:
			line=sam_file.readline()
			if not line:
				break
			else:
				columns = line.split("\t")
				gene = columns[2]
				gene = gene.strip()
				gene = gene.upper()
				if gene != '*':
					bowtie_al+=1
					umi = umi_list[sam_line_ind]
					barcode = corrected_bcs[sam_line_ind]
					if barcode:
						umi_bc_pair = (umi, barcode)
						if gene not in gene_to_umi_bc_dict[umi_bc_pair]:
							AS_score = int(columns[11][5:])
							AS_score_str = str(AS_score)
							dict_quality_scores[AS_score_str]+=1
							#if AS_score>=-3:
							sam_gene+=1
							gene_to_umi_bc_dict[umi_bc_pair].update({gene : 1})
							dict_quality[barcode]['high']+=1
							dict_barcode_counter[barcode] +=1
							barcode_counter+=1
							dict_genes_barcode[gene][barcode] += 1
							#else:
							#	dict_quality[barcode]['low']+=1
							#	sam_star += 1
						else:
							gene_to_umi_bc_dict[umi_bc_pair][gene] += 1
							dis_redund += 1
					else:
						read_has_no_bc +=1
				else:
					sam_star += 1
			sam_line_ind += 1

		# compute read count and alignment stats	
		alignment_score=sam_gene/(sam_gene+sam_star)
		alignment_score=round(alignment_score*100)
		bowtie_score=bowtie_al/(bowtie_al+sam_star)
		bowtie_score=round(bowtie_score*100)
		sam_file.close()
		print "Data stored in dictionaries........................................"
		print "Creating genes-cells matrix...\n"
		print gene_counter, "genes"
		print barcode_counter-1, "counted reads"
	return [gene_to_umi_bc_dict, dict_quality, dict_barcode_counter, 
			dict_genes_barcode, dis_redund, read_has_no_bc, sam_star, 
			bowtie_score, dict_quality_scores, sam_line_ind]	

############################### PREPROCESSING ###############################
print "\n"
print "***********************************"
print "** Loading files & cell barcodes **"
print "***********************************"
print "\n"

# Get all predefined groups: 

#Get a list of unprocessed fastq.gz files from the directory
pckl_files = [filename for filename in glob.glob(dir_path_pickle+ '*.pckl')]
#Sorting the list alphabetically in order to get R1 and R2 pairs) 
pckl_files.sort()

allgroups=list()
for p in pckl_files:
	p_file = open(p,'r')
	for i in range(0, 5):
		curr_groups = pickle.load(p_file)
	allgroups.extend(curr_groups)

# make sure that the barcodes don't overlap
top_bcs = [x[0] for x in allgroups]

D = defaultdict(list)
for i,item in enumerate(top_bcs):
  D[item].append(i)
D = {k:v for k,v in D.items() if len(v)>1}

for d in D.keys():
	for j in range(0, len(D[d])):
		if j==0: 
			numstring =''
		else:
			numstring =str(j)
		newbc = d + numstring
		allgroups[D[d][j]][0]=newbc

top_bcs = [x[0] for x in allgroups]

curr_fastq = 2

dir_path_alignment = dir_path_barcodes + 'alignments/'
dir_path_fastqs = dir_path_barcodes + 'fastqs/'

#### Name of data matrix to get ordering of barcodes from 
# matrix_files = [filename for filename in glob.glob(dir_path_barcodes +'*_matrix.txt') if 'pBC' not in filename]
# datafile = matrix_files[0]

#Get a list of unprocessed fastq.gz files from the directory
fastq_files = [filename for filename in glob.glob(dir_path_fastqs +'*.fastq.gz') if 'processed' not in filename]
#Sorting the list alphabetically in order to get R1 and R2 pairs) 
fastq_files.sort()
nb_fastqs = len(fastq_files)

for f1, f2 in grouped(fastq_files, 2):

	umi_list = []
	pre_barcode_list = []
	f1_base_name = '/'.join(f1.split('/')[:-1]) + '/' + f1.split(os.extsep)[0].split('/')[-1].split('_')[0]
	name = datafile.split('/')[-1].split('.')[0].split('_')[0]

	preprocess_results = preprocess_barcodes(f1, f2)
	umi_list = preprocess_results[0]
	pre_barcode_list = preprocess_results[1]
	total_reads = preprocess_results[2]

	#################################################################
	##### CELL BARCODE CORRECTION ###################################
	#################################################################

	# iterate over all barcodes in group and assemble pre_to_post_bc dictionary
	pre_2_post_bc = defaultdict(str)
	for i in range(0, len(allgroups)):
		topBC = allgroups[i][0]
		num_in_group = len(allgroups[i])
		for k in range(0, num_in_group):
			curr_bc = allgroups[i][k]
			pre_2_post_bc[curr_bc] = topBC

	# iterate over all the original barcodes and correct them in the list
	print "Assigning all barcodes to closest barcode"
	post_barcode_list = []
	uncorrected = 0
	corrected = 0
	corrected_bcs_file=open(f1_base_name+'_barcode.txt', 'w+', 1)	
	for pre_barcode in pre_barcode_list:
		corrected_bc = pre_2_post_bc[pre_barcode] # will be '' if its not in the list
		corrected_bcs_file.write(corrected_bc+'\n')
		post_barcode_list.append(corrected_bc)
		if (pre_barcode != corrected_bc) and corrected_bc:
			corrected +=1
		if corrected_bc == '':
			uncorrected += 1
	corrected_bcs_file.close()

	print '\t', total_reads, 'Total barcode reads'
	print "..................................................................."
	os.system("gzip "+f1_base_name+'_processed.fastq')

	# read fasta to get gene names
	fasta_return_dicts = get_gene_names_from_fasta(fasta_path)
	dict_gene_names = fasta_return_dicts[0]
	dict_gene_counter = fasta_return_dicts[1]
	gene_counter = len(dict_gene_counter)

	start_time = time.time()
	#Listing *_processed.fastq files:
	preprocessed_files = os.listdir(dir_path_barcodes)
	#Cleaning list of files in the directory (Only *processed.fastq files):
	files_processed = [f for f in preprocessed_files if 'processed.fastq.gz' in f]
	files_processed.sort()
	bowtie_opt = ' '.join(bowtie2_options)

	######################### align reads with bowtie2 ##########################
	print "\n"
	print "**********************************"
	print "**      Starting alignment      **"
	print "**********************************"
	print "\n"
	
	sam_file_name = align_reads(f1_base_name+'_processed.fastq.gz')

	read_sam_list = read_sam_umi(sam_file_name, f1_base_name)		
	gene_to_umi_bc_dict =read_sam_list[0]
	dict_quality = read_sam_list[1]
	dict_barcode_counter = read_sam_list[2]
	dict_genes_barcode = read_sam_list[3]
	dis_redund = read_sam_list[4]
	read_has_no_bc = read_sam_list[5]
	sam_star = read_sam_list[6]
	bowtie_score = read_sam_list[7]
	dict_quality_scores = read_sam_list[8]
	total_reads = read_sam_list[9]

	# 	Settle the UMIs
	for key in gene_to_umi_bc_dict.keys():
		entry = gene_to_umi_bc_dict[key] # also a dictionary
		sortedpBCs = sorted(entry)
		vals = [entry[x] for x in sortedpBCs]
		indexes = [i for i, x in enumerate(vals) if x==max(vals)]
		if (len(indexes)>1):
			gene_to_umi_bc_dict[key]='none'
		else:
			maxgene = sortedpBCs[indexes[0]]
			gene_to_umi_bc_dict[key]=maxgene

	cellbc_to_gene_dict = defaultdict(lambda: defaultdict(int))

	# Assemble cell barcode - gene count matrix:
	for barcode in top_bcs:
		curr_keys = [pair for pair in gene_to_umi_bc_dict.keys() if pair[1]==barcode]
		curr_keys.sort()
		genevals = [gene_to_umi_bc_dict[key] for key in curr_keys]
		for key_gene in dict_gene_names:
			gene_count=genevals.count(key_gene)
			cellbc_to_gene_dict[barcode].update({key_gene:gene_count})

	d = open(datafile, 'r')
	fin_bcs = list(d.readline().split('\n')[0].split('\t')[1:])
	bc_2_column = {bc:i+1 for (i,bc) in enumerate(fin_bcs)}

	# remove barcodes from cellbc_to_gene_dict that aren't in the final matrix file:
	zerobcs = list(set(cellbc_to_gene_dict.keys()) - set(fin_bcs))

	for zerobc in zerobcs:
		del cellbc_to_gene_dict[zerobc]

	matrix = [[0 for x in range(len(fin_bcs)+1)] for y in range(gene_counter+1)]

	# write barcode header for each column
	for bc in fin_bcs:
		col_num = bc_2_column[bc]
		matrix[0][col_num] = bc

	# write the gene name header for each row
	for key_gene in dict_gene_names:
		row_num = dict_gene_counter[dict_gene_names[key_gene]]
		matrix[row_num][0] = dict_gene_names[key_gene]

	# write read counts for each gene at each barcode
	reads_counted = 0
	for key_barcode in cellbc_to_gene_dict.keys():
		col_num = bc_2_column[key_barcode]
		for key_gene in cellbc_to_gene_dict[key_barcode].keys():
			row_num = dict_gene_counter[dict_gene_names[key_gene]]
			matrix[row_num][col_num] += cellbc_to_gene_dict[key_barcode][key_gene]
			reads_counted += cellbc_to_gene_dict[key_barcode][key_gene]

	print "Genes-cells matrix created........................................."

	matrix_file = open(dir_path_barcodes + name + '_pBC_matrix.txt', 'w+')
	for item in matrix:
			matrix_file.write('\t'.join([str(i) for i in item])+'\n')
	matrix_file.close()

# def getKey(item):
# 	return(item[1])

# sortedpairs = sorted(gene_to_umi_bc_dict.keys(), key=getKey)