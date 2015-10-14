# -*-coding:Utf-8 -*

#######################################
## plate.py is designed to preprocess
## the fastq files from a plate exper 
## iment, and align data to a referen
## ce genome.
#######################################
## Developed by:
## Paul Rivaud
## paulrivaud.info@gmail.com
## Summer 2015
#######################################

from __future__ import division
from itertools import izip
from collections import defaultdict
from plate_params import *
import re
import math
import os
import os.path
import time
import gzip

def round_figures(x, n):
	return round(x, int(n - math.ceil(math.log10(abs(x)))))

'''
#Getting a files list from the directory
files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory(removing *noTA.fastq and *umi.txt files):
files_noNoTA = [f for f in files if 'noTA' not in f]
files_noNoTA_noUMI = [f for f in files_noNoTA if 'umi' not in f]
#Removing files that do not have a .fastq file extension:
fastq_files = [f for f in files_noNoTA_noUMI if '.fastq.gz' in f]
#Sorting the list alphabetically in order to get R1 and R2 pairs) 
fastq_files.sort()
############################### PREPROCESSING ###############################
print "\n"
print "**********************************"
print "**    Starting preprocessing    **"
print "**********************************"
print "\n"

nb_fastqs = len(fastq_files)
curr_fastq = 1

for f1 in fastq_files:
	percent = int((curr_fastq/nb_fastqs)*100)
	curr_fastq+=1
	if os.path.isfile(dir_path_fastqs+f1.split(os.extsep)[0]+'_noTA.fastq.gz'):
		print f1
		print "\tFiles already preprocessed. Moving on to other files...",percent,"%"
	else:
		modified = True
		print "Begin..."
		print "Opening fastq file..."
		fastq1_path = dir_path_fastqs+f1.split(os.extsep)[0]
		file1_fastq = gzip.open(fastq1_path+'.fastq.gz','rb')
		#Dictionary containing (seq)-(umis list) pairs
		seq_dictionary = defaultdict(list)
		print f1
		print "\tReading file..."
		print "\tWriting file..."
		file_noTA = gzip.open(fastq1_path+'_noTA.fastq.gz', 'wb')
		#Stats about the trimming process
		total_reads = 0
		while True:
			f1_line1 = file1_fastq.readline()
			f1_line2 = file1_fastq.readline()
			f1_line3 = file1_fastq.readline()
			f1_line4 = file1_fastq.readline()
			if not f1_line1:
				break
			else:
				total_reads+=1
			if tso not in f1_line2:
				#Removing polyTs
				index_t =  [m.start() for m in re.finditer('(?='+error_string_t+')', f1_line2)]
				reversed_index_t = list(reversed(index_t))
				for x in reversed_index_t:
					while f1_line2[x] == 'T':
						f1_line2 = f1_line2[:x] + f1_line2[x+1:]
						f1_line4 = f1_line4[:x] + f1_line4[x+1:]
				#Removing polyAs
				index_a =  [m.start() for m in re.finditer('(?='+error_string_a+')', f1_line2)]
				reversed_index_a = list(reversed(index_a))
				for x in reversed_index_a:
					while f1_line2[x] == 'A':
						f1_line2 = f1_line2[:x] + f1_line2[x+1:]
						f1_line4 = f1_line4[:x] + f1_line4[x+1:]
				#Checking trimmed sequence length
				if len(f1_line2) >= minimum_length:
					file_noTA.write(f1_line1)
					file_noTA.write(f1_line2)
					file_noTA.write(f1_line3)
					file_noTA.write(f1_line4)
		file1_fastq.close()
		file_noTA.close()
		print "\tTrimming results:"
		print "\tTotal reads: ",total_reads
		print "...................................................................",percent,"%"

'''
############################### ALIGNMENT ###############################
start_time = time.time()
#Listing *_noTA.fastq files:
preprocessed_files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory (Only *noTA.fastq files):
files_noTA = [f for f in preprocessed_files if 'fastq.gz' in f]
files_noTA.sort()
bowtie_opt = ' '.join(bowtie2_options)
print "\n"
print "**********************************"
print "**      Starting alignment      **"
print "**********************************"
print "\n"

nb_noTA = len(files_noTA)
curr_noTA = 1
print files_noTA
for file_noTA in files_noTA:
	percent = int((curr_noTA/nb_noTA)*100)
	curr_noTA+=1
	if os.path.isfile(dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam'):
		print "Skipping alignment: "+file_noTA.split(os.extsep)[0]+'.sam already exists...',percent,"%"
	else:
		modified = True
		print "Starting alignment with bowtie2..."
		os.system("nice "+bowtie2_dir+" bowtie2 "+bowtie_opt+\
			" -x "+reference_genome+" -U "+dir_path_fastqs+file_noTA+\
			" -S "+dir_path_alignment+file_noTA.split(os.extsep)[0]+'.sam')
		print file_noTA+" aligned..."
		print file_noTA.split(os.extsep)[0]+'.sam created...'
		print "...................................................................",percent,"%"
total_time = time.time() - start_time
print "Reads alignment time:"
print int(total_time/60),"min",int(total_time%60),"sec"

print "\n"
print "**********************************"
print "**       Retrieving data        **"
print "**********************************"
print "\n"

sum_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])
if os.path.isfile(sum_path+'matrix.txt') and modified == False:
	print "matrix.txt already exists... 100%"
else:
	#Listing .sam files:
	sam_files = os.listdir(dir_path_alignment)
	sam_files.sort()
	nb_sam = len(sam_files)
	curr_sam = 1
	gene_counter = 1
	experiment_counter = 1
	dict_genes_exp = defaultdict(dict)
	dict_gene_names = defaultdict(dict)
	dict_gene_counter = defaultdict(int)
	dict_exp_counter = defaultdict(str)

	fasta_file = open(fasta_path,'r')
	while True:
		line = fasta_file.readline()
		if not line:
			break
		else:
			if line[:1] == '>':
				gene_symbol=line.split(' ')[1]
				gene_symbol=gene_symbol[5:]
				gene_symbol=gene_symbol.replace('\n','')
				if gene_symbol not in dict_gene_counter:
					dict_gene_counter[gene_symbol] = gene_counter
					gene_counter+=1
				gene_nm=line.split(' ')[0]
				gene_nm=gene_nm[1:]
				gene_nm=gene_nm.replace('\n','')
				dict_gene_names[gene_nm] = gene_symbol
	fasta_file.close()


	print "Checking SAM files..."
	for sam_file in sam_files:
		current_sam = open(dir_path_alignment+sam_file,'r')
		print"Opened:",sam_file
		dict_exp_counter[experiment_counter] = sam_file
		while True:
			line=current_sam.readline()
			if not line:
				break
			else:
				columns = line.split("\t")
				gene = columns[2]
				#If read aligned, columns[2] is different from '*'
				if gene != '*':
					if gene in dict_genes_exp:
						if experiment_counter in dict_genes_exp[gene].keys():
							dict_genes_exp[gene][experiment_counter] +=1
						else:
							dict_genes_exp[gene][experiment_counter] = 1
					else:
						dict_genes_exp[gene] = {experiment_counter : 1}

		experiment_counter+=1
		current_sam.close()
		percent = int((curr_sam/nb_sam)*100)
		curr_sam+=1
		print "...................................................................",percent,"%"
	print "Data stored in dictionaries........................................",percent,"%"
	print "Creating genes-cells matrix...\n"
	print gene_counter, "genes"
	print experiment_counter-1, "experiments"

	matrix = [[0 for x in range(experiment_counter)] for x in range(gene_counter)]
	#Fill col names:
	for key_experiment in dict_exp_counter:
		matrix[0][key_experiment] = dict_exp_counter[key_experiment]
	for key_gene in dict_gene_names:
		row_num = dict_gene_counter[dict_gene_names[key_gene]]
		matrix[row_num][0] = dict_gene_names[key_gene]
	#Fill matrix values:
	for key_gene in dict_gene_names:
		row_num = dict_gene_counter[dict_gene_names[key_gene]]
		for key_experiment in dict_genes_exp[key_gene]:
			matrix[row_num][key_experiment]+=dict_genes_exp[key_gene][key_experiment]
	print "Genes-cells matrix created.........................................",percent,"%"
	matrix_file = open(sum_path+'matrix.txt', 'w+')
	for item in matrix:
		matrix_file.write('\t'.join([str(i) for i in item])+'\n')
	matrix_file.close()
	print "Done..."	
print "\n"
print "**********************************"
print "***********Pipeline end***********"
print "**********************************"
print "\n"
