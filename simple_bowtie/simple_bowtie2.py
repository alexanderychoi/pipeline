# -*-coding:Utf-8 -*
from simple_bowtie2_params import *
import gzip
import os
import os.path
import time
import re

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

for f in fastq_files:
	percent = int((curr_fastq/nb_fastqs)*100)
	curr_fastq+=1
	if os.path.isfile(dir_path_fastqs+f.split(os.extsep)[0]+'_noTA.fastq.gz'):
		print f
		print "\tFile already preprocessed. Moving on to other file...",percent,"%"
	else:
		print "Begin..."
		print "Opening fastq file..."
		fastq1_path = dir_path_fastqs+f.split(os.extsep)[0]
		file1_fastq = gzip.open(fastq1_path+'.fastq.gz','rb')
		#Stats about the trimming process
		total_reads = 0
		#Dictionary containing (seq)-(umis list) pairs
		print f
		print "\tReading file..."
		print "\tWriting file..."
		file_noTA = gzip.open(fastq1_path+'_noTA.fastq.gz', 'wb')
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
		print "\tTotal reads: ",total_reads
		print "...................................................................",percent,"%"

############################### ALIGNMENT ###############################
start_time = time.time()
#Listing *_noTA.fastq files:
preprocessed_files = os.listdir(dir_path_fastqs)
#Cleaning list of files in the directory (Only *noTA.fastq files):
files_noTA = [f for f in preprocessed_files if 'noTA.fastq.gz' in f]
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
	if os.path.isfile(dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam.gz'):
		print "Skipping alignment: "+os.path.splitext(file_noTA)[0]+'.sam.gz already exists...',percent,"%"
	else:
		print "Starting alignment with bowtie2 for", file_noTA
		os.system("nice "+bowtie2_dir+" bowtie2 "+bowtie_opt+\
			" -x "+reference_genome+" -U "+dir_path_fastqs+file_noTA+\
			" -S "+dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam')
		print file_noTA+" aligned..."
		print os.path.splitext(file_noTA)[0]+'.sam created...'
		print "...................................................................",percent,"%"
		os.system("gzip "+dir_path_alignment+os.path.splitext(file_noTA)[0]+'.sam')
		print "Sam file compressed..."
		print "...................................................................",percent,"%"
total_time = time.time() - start_time
print "Reads alignment time:"
print int(total_time/60),"min",int(total_time%60),"sec"


print "\n"
print "**********************************"
print "***********Pipeline end***********"
print "**********************************"
print "\n"