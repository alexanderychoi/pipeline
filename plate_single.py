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

def round_figures(x, n):
	return round(x, int(n - math.ceil(math.log10(abs(x)))))

modified = False

#Getting a files list from the directory
files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory(removing *noTA.fastq and *umi.txt files):
files_noNoTA = [f for f in files if 'noTA' not in f]
files_noNoTA_noUMI = [f for f in files_noNoTA if 'umi' not in f]
#Removing files that do not have a .fastq file extension:
fastq_files = [f for f in files_noNoTA_noUMI if '.fastq' in f]
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
	if os.path.isfile(dir_path_fastqs+os.path.splitext(f1)[0]+'_noTA.fastq'):
		print f1
		print "\tFiles already preprocessed. Moving on to other files...",percent,"%"
	else:
		modified = True
		print "Begin..."
		print "Opening fastq file..."
		fastq1_path = dir_path_fastqs+os.path.splitext(f1)[0]
		file1_fastq = open(fastq1_path+'.fastq','r')
		#Dictionary containing (seq)-(umis list) pairs
		seq_dictionary = defaultdict(list)
		print f1
		print "\tReading file..."
		print "\tWriting file..."
		file_noTA = open(fastq1_path+'_noTA.fastq', 'w+')
		#Stats about the trimming process
		total_reads = 0
		saved_reads = 0
		dismissed_reads = 0
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
				else:
					dismissed_reads+=1
			else:
				dismissed_reads+=1
		file1_fastq.close()
		file2_fastq.close()
		file_noTA.close()
		print "\tTrimming results:"
		print "\tTotal reads: ",total_reads
		print "\tReads saved: ", round_figures((saved_reads*100)/total_reads,3), "%"
		print "\tReads dismissed: ", round_figures((dismissed_reads*100)/total_reads,3), "%"
		print "...................................................................",percent,"%"


############################### ALIGNMENT ###############################
start_time = time.time()
#Listing *_noTA.fastq files:
preprocessed_files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory (Only *noTA.fastq files):
files_noTA = [f for f in preprocessed_files if 'noTA' in f]
files_noTA.sort()
bowtie_opt = ' '.join(bowtie2_options)
print "\n"
print "**********************************"
print "**      Starting alignment      **"
print "**********************************"
print "\n"

nb_noTA = len(files_noTA)
curr_noTA = 1

for file_noTA in files_noTA:
	percent = int((curr_noTA/nb_noTA)*100)
	curr_noTA+=1
	if os.path.isfile(dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam'):
		print "Skipping alignment: "+os.path.splitext(file_noTA)[0]+'.sam already exists...',percent,"%"
	else:
		modified = True
		print "Starting alignment with bowtie2..."
		os.system("nice "+bowtie2_dir+" bowtie2 "+bowtie_opt+\
			" -x "+reference_genome+" -U "+dir_path_fastqs+file_noTA+\
			" -S "+dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam')
		print file_noTA+" aligned..."
		print os.path.splitext(file_noTA)[0]+'.sam created...'
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
if os.path.isfile(sum_path+'genes.txt') and modified == False:
	print "genes.txt already exists... 100%"
else:
	#Listing .sam files:
	sam_files = os.listdir(dir_path_alignment)
	sam_files.sort()
	nb_sam = len(sam_files)
	curr_sam = 1
	file_genes = open(sum_path+'genes.txt', 'w+')
	gene_counter = 0
	experiment_counter = 1
	dict_barcode_genes = defaultdict(dict)
	dict_gene_counter = defaultdict(int)

	print "Checking SAM files..."
	for sam_file in sam_files:
		current_sam = open(dir_path_alignment+sam_file,'r')
		while True:
			line=sam_file.readline()
			if not line:
				break
			else:
				columns = line.split("\t")
				gene = columns[2]
				#If read aligned, columns[2] is different from '*'
				if gene != '*':
					if gene not in dict_gene_counter:
						dict_gene_counter[gene] = gene_counter
						gene_counter+=1
					if gene in dict_barcode_genes[experiment_counter].keys():
							dict_barcode_genes[experiment_counter][gene] +=1
						else:
							dict_barcode_genes[experiment_counter][gene] = 1
					else:
						dict_barcode_genes[experiment_counter] = {gene : 1}
		experiment_counter+=1

		for line in current_sam.readlines():
			columns = line.split("\t")
			gene = columns[2]
			if gene in gene_dictionary:
				gene_dictionary[gene] += 1
			else:
				gene_dictionary[gene] = 1
		current_sam.close()
		percent = int((curr_sam/nb_sam)*100)
		curr_sam+=1
		print percent,"%"
	for key in gene_dictionary:
		file_genes.write(key+"\t"+str(gene_dictionary[key])+"\n")
	file_genes.close()
	print "Done..."
print "Pipeline end reached..."
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	while True:
		line=sam_file.readline()
		if not line:
			break
		else:
			columns = line.split("\t")
			gene = columns[2]
			barcode = barcode_file.readline()
			barcode = barcode.replace('\n','')
			#If read aligned, columns[2] is different from '*'
			if gene != '*':
				if gene not in dict_gene_counter:
					dict_gene_counter[gene] = gene_counter
					gene_counter+=1
				if barcode in dict_barcode_genes:
					if gene in dict_barcode_genes[barcode].keys():
						dict_barcode_genes[barcode][gene] +=1
					else:
						dict_barcode_genes[barcode][gene] = 1
				else:
					dict_barcode_genes[barcode] = {gene : 1}
					barcode_counter+=1
	barcode_file.close()
	print "Data stored in dictionaries........................................",percent,"%"
	print "Creating genes-cells matrix...\n"

	print gene_counter, "genes"
	print barcode_counter, "cells"
	matrix = [[0 for x in range(barcode_counter+1)] for x in range(gene_counter+1)]
	col_num=1
	for key_barcode in dict_barcode_genes:
		matrix[0][col_num] = key_barcode
		for key_gene in dict_barcode_genes[key_barcode]:
			col_gene = dict_gene_counter[key_gene]+1
			matrix[col_gene][0] = key_gene
			matrix[col_gene][col_num]=dict_barcode_genes[key_barcode][key_gene]
		col_num+=1
	print "Genes-cells matrix created.........................................",percent,"%"
	matrix_file = open(sum_path+'matrix.txt', 'w+')
	for item in matrix:
		matrix_file.write('\t'.join([str(i) for i in item])+'\n')
	matrix_file.close()
print "Ready for PCA..."
print "\n"
print "**********************************"
print "***********Pipeline end***********"
print "**********************************"
print "\n"
