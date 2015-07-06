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
from plate_params_kallisto import *
import re
import math
import os
import os.path
import time

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def round_figures(x, n):
	return round(x, int(n - math.ceil(math.log10(abs(x)))))

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

avg_lg = 0


for f1, f2 in grouped(fastq_files, 2):
	if os.path.isfile(dir_path_fastqs+os.path.splitext(f1)[0]+'_noTA.fastq'):
		print f1
		print f2
		print "\tFiles already preprocessed. Moving on to other files..."
	else:
		pass
		print "Begin..."
		print "Opening fastq files..."
		fastq1_path = dir_path_fastqs+os.path.splitext(f1)[0]
		fastq2_path = dir_path_fastqs+os.path.splitext(f2)[0]
		file1_fastq = open(fastq1_path+'.fastq','r')
		file2_fastq = open(fastq2_path+'.fastq','r')
		#Dictionary containing (seq)-(umis list) pairs
		seq_dictionary = defaultdict(list)

		print f1
		print f2
		print "\tReading files..."
		print "\tWriting files..."
		file_noTA = open(fastq1_path+'_noTA.fastq', 'w+')
		file_umi = open(fastq1_path+'_umi.txt', 'w+')
		#Stats about the trimming process
		total_reads = 0
		saved_reads = 0
		dismissed_reads = 0
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

		while True:
			f1_line1 = file1_fastq.readline()
			f1_line2 = file1_fastq.readline()
			f1_line3 = file1_fastq.readline()
			f1_line4 = file1_fastq.readline()
			f2_line1 = file2_fastq.readline()
			f2_line2 = "".join(complement.get(base, base) for base in reversed(file2_fastq.readline()))
			f2_line3 = file2_fastq.readline()
			f2_line4 = "".join(reversed(file2_fastq.readline()))
			if not f1_line1:
				break
			else:
				total_reads+=1

			if tso not in f2_line2 and tso not in f2_line2:
				if f1_line2[0:tac_length] == 'TAC':
					umi = f1_line2[tac_length:tac_length+umi_length]
					if umi != 'TTTTT':
						cat_seq = f1_line2[tac_length+umi_length+g_length:-1]+f2_line2[1:]+'\n'
						cat_qual = f1_line4[tac_length+umi_length+g_length:-1]+f2_line4[1:]+'\n'
						#Removing polyTs
						index_t =  [m.start() for m in re.finditer('(?='+error_string_t+')', cat_seq)]
						reversed_index_t = list(reversed(index_t))
						for x in reversed_index_t:
							while cat_seq[x] == 'T':
								cat_seq = cat_seq[:x] + cat_seq[x+1:]
								cat_qual = cat_qual[:x] + cat_qual[x+1:]
						#Removing polyAs
						index_a =  [m.start() for m in re.finditer('(?='+error_string_a+')', cat_seq)]
						reversed_index_a = list(reversed(index_a))
						for x in reversed_index_a:
							while cat_seq[x] == 'A':
								cat_seq = cat_seq[:x] + cat_seq[x+1:]
								cat_qual = cat_qual[:x] + cat_qual[x+1:]
						#Checking trimmed sequence length
						if len(cat_seq) >= minimum_length:
							if cat_seq in seq_dictionary:
								if umi not in seq_dictionary[cat_seq]:
									seq_dictionary[cat_seq].append(umi)
									file_umi.write(umi+'\n')
									file_noTA.write(f1_line1)
									file_noTA.write(cat_seq)
									file_noTA.write(f1_line3)
									file_noTA.write(cat_qual)
								else:
									dismissed_reads+=1
							else:
								seq_dictionary[cat_seq].append(umi)
								file_umi.write(umi+'\n')
								file_noTA.write(f1_line1)
								file_noTA.write(cat_seq)
								file_noTA.write(f1_line3)
								file_noTA.write(cat_qual)
								saved_reads+=1
						else:
							dismissed_reads+=1
					else:
						dismissed_reads+=1
				else:
					dismissed_reads+=1
			else:
				dismissed_reads+=1

		file_umi.close()
		file1_fastq.close()
		file2_fastq.close()
		file_noTA.close()

		print "\tTrimming results:"
		print "\tTotal reads: ",total_reads
		print "\tReads saved: ", round_figures((saved_reads*100)/total_reads,3), "%"
		print "\tReads dismissed: ", round_figures((dismissed_reads*100)/total_reads,3), "%"


############################### ALIGNMENT ###############################
#Listing *_noTA.fastq files:
preprocessed_files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory (Only *noTA.fastq files):
files_noTA = [f for f in preprocessed_files if 'noTA' in f]
files_noTA.sort()
print "\n"
print "**********************************"
print "**      Starting alignment      **"
print "**********************************"
print "\n"
print "Checking for index..."
if os.path.isfile(reference_idx):
	print "Index already exists..."
else:
	print "Index not found. Creating index..."
	os.system("kallisto index -i "+reference_idx+" "+reference_genome)
	print "Index created..."

for file_noTA in files_noTA:
	if os.path.isfile(dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam'):
		print "Skipping alignment: "+os.path.splitext(file_noTA)[0]+'.sam already exists.'
	else:
		print "Starting alignment with Kallisto..."
		os.system("kallisto quant -i "+reference_idx+" -o out --pseudobam --single -l 140 "+dir_path_fastqs+file_noTA+" > "+dir_path_alignment+os.path.splitext(file_noTA)[0]+".sam")
		print file_noTA+" aligned..."
		print os.path.splitext(file_noTA)[0]+'.sam created...'
