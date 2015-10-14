# -*-coding:Utf-8 -*
## TITAN WILL LIVE FOR EVER IN OUR HEARTS
#######################################
## drop.py is designed to preprocess
## the fastq files from a plate exper 
## iment, and align data to a referen
## ce genome.
#######################################
## Developed by:
## Paul Rivaud
## paulrivaud.info@gmail.com
## Summer 2015
## Maintained by:
## Graham Heimberg
## Graham.Heimberg@ucsf.edu
#######################################

from __future__ import division
from itertools import izip
from collections import defaultdict
from drop_titan_params import *
from generate_html_report import *
import re
import math
import os
import os.path
import time
import gzip
import operator
import glob
import numpy as np
from matplotlib import pyplot as plt
from scipy import cluster, io
import sys

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

def longest_common_substring(s1, s2):
	# used for finding large chunks of the TSO 
	# oligo sequence within the read
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return x_longest, longest

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def str_dist(s1, s2):
	# dynamic programming algorithm to compute the 
	# edit distance (levenshtein distance) between two strings
	if s1==s2:
		return 0
	previous_row = range(len(s2) + 1)
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
			deletions = current_row[j] + 1       # than s2
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]

def compute_bc_dist_mat(bc_list):
	# Compute a distance matrix using the levenshtein distance
	# metric from a list of barcodes
	bc_num = len(bc_list)
	bc_len = len(bc_list[0])
	dist_mat = np.zeros([bc_num,bc_num])
	for i in range(len(bc_list)):
		for j in range(i, len(bc_list)):
			dist_mat[i][j] = str_dist(bc_list[i], bc_list[j])
	dist_mat = dist_mat + dist_mat.transpose()
	#dist_mat[dist_mat>5] = bc_len
	return dist_mat

def measure_dist(bc, centroids_list):
	# compute distance of a barcode to each centroid in a list
	dist_to_centroids=[0 for i in range(len(centroids_list))]
	for i in range(len(centroids_list)):
		dist_to_centroids[i]=str_dist(bc,centroids_list[i])
	return dist_to_centroids

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

def preprocess_fastq_file_pair(f1, f2):
	print "processing fastq file pairs"
	print "Opening fastq files..."
	file1_fastq = gzip.open(f1,'rb')
	file2_fastq = gzip.open(f2,'rb')
	print f1
	print f2
	print "\tReading files..."
	print "\tWriting files..."
	file_processed = open(f1_base_name+'_processed.fastq', 'w+', 1000)
	file_umi = open(f1_base_name+'_umi.txt', 'w+',1)
	file_pre_barcode = open(f1_base_name+'_pre_barcode.txt', 'w+', 1000)
	#Stats about the trimming process
	total_reads = 0
	preprocessing_saved_reads = 0
	dismissed_reads = 0
	dis_tso = 0
	dis_no_tac = 0
	tac_length = len(str_search)
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
		if (str_search in f1_line2[:tac_length]) and f1_line2[:(tac_length+3)] != (str_search + 'GGG'):
			if tso not in f2_line2:
				curr_read = fastq_read(f2_line1, f2_line2, f2_line3, f2_line4, f1_line2)
				curr_read.remove_polyA()
				#curr_read.remove_tso()
				if curr_read.line1:
					file_umi.write(curr_read.umi+'\n')
					umi_list.append(curr_read.umi)
					pre_barcode_list.append(curr_read.pre_barcode)					
					barcode_count_dict[curr_read.pre_barcode] += 1
					file_pre_barcode.write(curr_read.pre_barcode+'\n')				
					curr_read.write_read(file_processed)
					preprocessing_saved_reads+=1
				else:
					dismissed_reads+=1					
			else:
				dis_tso+=1
				dismissed_reads+=1
		else:
			dis_no_tac+=1
			dismissed_reads+=1
	file_pre_barcode.close()
	file_umi.close()		
	file1_fastq.close()
	file2_fastq.close()
	file_processed.close()
	return [umi_list, pre_barcode_list, barcode_count_dict, preprocessing_saved_reads, dismissed_reads, dis_tso, dis_no_tac, total_reads]

def align_reads(fastq_file):
	start_time = time.time()
	fastq_file_name = fastq_file.split('/')[-1]
	alignment_file_name = dir_path_alignment + fastq_file.split('/')[-1].split(os.extsep)[0] + '.sam'
	if os.path.isfile(fastq_file.split(os.extsep)[0]+'.sam.gz'):
		print "Skipping alignment: " + fastq_file.split(os.extsep)[0] + '.sam.gz already exists...'
	else:
		print "Starting alignment with bowtie2 for", fastq_file_name
		os.system("nice " + bowtie2_dir + " bowtie2 " + bowtie_opt +\
			" -x " + reference_genome + " -U " + fastq_file +\
			" -S " + alignment_file_name)
		print 
		os.system("gzip " + alignment_file_name)
		print fastq_file_name + "aligned and file compressed"
		print "..................................................................."
	total_time = time.time() - start_time
	print "Reads alignment time:"
	print int(total_time/60),"min",int(total_time%60),"sec"
	return alignment_file_name

def txt_file_to_list(filename):
	in_file = open(filename, 'r', 1000)
	row_list = in_file.read().splitlines()
	inf_file.close()
	return row_list

def read_sam(sam_file_name):
	############################### Data gathering ###############################
	print "\n"
	print "**********************************"
	print "**       Retrieving data        **"
	print "**********************************"
	print "\n"

	if os.path.isfile(common_path+'matrix.txt'):
		print "matrix.txt already exists.......................................... 100 %"
	else:
		barcode_counter = 1
		sam_file = gzip.open(sam_file_name + '.gz', 'rb')
		list_f = os.listdir(dir_path_fastqs)
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
		gene_to_umi_bc_dict = defaultdict(lambda: set())

		while True:
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
					pre_barcode = post_barcode_list[sam_line_ind]
					barcode = pre_2_post_bc[pre_barcode]
					if barcode:
						umi_bc_pair = (umi, barcode)
						if umi_bc_pair not in gene_to_umi_bc_dict[gene]:
							AS_score = int(columns[11][5:])
							AS_score_str = str(AS_score)
							dict_quality_scores[AS_score_str]+=1
							#if AS_score>=-3:
							sam_gene+=1
							gene_to_umi_bc_dict[gene].add(umi_bc_pair)
							dict_quality[barcode]['high']+=1
							dict_barcode_counter[barcode] +=1
							barcode_counter+=1
							dict_genes_barcode[gene][barcode] += 1
							#else:
							#	dict_quality[barcode]['low']+=1
							#	sam_star += 1
						else:
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
		print gene_counter-1, "genes"
		print barcode_counter-1, "counted reads"
	return [gene_to_umi_bc_dict, dict_quality, dict_barcode_counter, 
			dict_genes_barcode, dis_redund, read_has_no_bc, sam_star, 
			bowtie_score, dict_quality_scores]	

def write_to_mat_file(matrix_list, filename):
	# save MATLAB variable with each barcode
	centroid_barcodes = np.zeros(cell_num, dtype=object)
	centroid_barcodes[:] = matrix[0][1:]
	# save MATLAB variable with each gene name
	gene_names = np.zeros(len(matrix)-1, dtype=object)
	gene_names[:] = [matrix[i][0] for i in range(1,len(matrix))]
	# save MATLAB variable with each 
	read_counts = np.array(matrix, dtype=object)
	read_counts = np.array(read_counts[1:, 1:], dtype=float)
	io.savemat(filename, mdict={'centroid_barcodes':centroid_barcodes, 
						'gene_names':gene_names, 'read_counts':read_counts})

#Get a list of unprocessed fastq.gz files from the directory
fastq_files = [filename for filename in glob.glob(dir_path_fastqs + '*.fastq.gz') if 'processed' not in filename]
#Sorting the list alphabetically in order to get R1 and R2 pairs) 
fastq_files.sort()
############################### PREPROCESSING ###############################
print "\n"
print "**********************************"
print "**    Starting preprocessing    **"
print "**********************************"
print "\n"

nb_fastqs = len(fastq_files)
curr_fastq = 2
umi_list = []
pre_barcode_list = []
barcode_count_dict = defaultdict(int)
common_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])

# initialize variables for saving to .mat files
matrix_list = []
matrix_name_list = []

for f1, f2 in grouped(fastq_files, 2):
	curr_fastq+=2
	f1_base_name = f1.split(os.extsep)[0]

	if os.path.isfile(dir_path_fastqs+f1.split(os.extsep)[0]+'_processed.fastq.gz'):
		print f1
		print f2
		pre_barcode_list = txt_file_to_list(f1_base_name + '_pre_barcode.txt')
		umi_list = txt_file_to_list(f1_base_name + '_umi.txt')
		print "\tFiles already preprocessed. Loading in parameters ..."

	else:
		preprocess_results = preprocess_fastq_file_pair(f1, f2)
		umi_list = preprocess_results[0]
		pre_barcode_list = preprocess_results[1]
		barcode_count_dict = preprocess_results[2]
		preprocessing_saved_reads = preprocess_results[3]
		dismissed_reads = preprocess_results[4]
		dis_tso = preprocess_results[5]
		dis_no_tac = preprocess_results[6]
		total_reads = preprocess_results[7]

		#################################################################
		##### BARCODE CORRECTION ########################################
		#################################################################

		num_unique_pre_bcs = len(barcode_count_dict)
	
		print 'starting clustering'

		# sort barcodes on occurances to cluster top 200
		top_bc_num = cell_num*10;
		sorted_bcs = sorted(barcode_count_dict.items(), key=operator.itemgetter(1))
		top_bcs={}
		for i in range(max(num_unique_pre_bcs-top_bc_num, 0), num_unique_pre_bcs):
			top_bcs[sorted_bcs[i][0]] = sorted_bcs[i][1]

		bc_dist_mat = compute_bc_dist_mat(top_bcs.keys())

		bc_linkage = cluster.hierarchy.linkage(bc_dist_mat, method='complete')
		clusters = cluster.hierarchy.fcluster(bc_linkage, cell_num, criterion='maxclust')
		num_clust = np.max(clusters)

		# find the bc that has the lowest distance to all other bc's within a cluster
		centroids = []
		for i in range(1,num_clust+1):
			cluster_dist_mat = bc_dist_mat[clusters==i,:]
			cluster_dist_mat = cluster_dist_mat[:,clusters==i]
			cluster_bcs = [bc for j,bc in enumerate(top_bcs) if clusters[j]==i]
			centroids.append(cluster_bcs[np.argmin(np.sum(cluster_dist_mat, axis=0))])

		# iterate over unique barcodes to populate correction dictionary
		pre_2_post_bc = defaultdict(str)
		unique_pre_bcs = set(pre_barcode_list)
		recovered_bc = 0
		unrecovered_bc_list = []
		for bc in unique_pre_bcs:
			bc_dist = measure_dist(bc, centroids)
			min_ind = np.argmin(bc_dist)
			if bc_dist[min_ind]<3:
				pre_2_post_bc[bc] = centroids[min_ind]
				recovered_bc += 1
			else:
				unrecovered_bc_list.append(bc)

		# iterate over all the original barcodes and write all the corrected barcodes to a file
		print "Assigning all barcodes to closest barcode"
		corrected_bcs_file=open(f1_base_name+'_barcode.txt', 'w+', 1)		
		post_barcode_list = []
		unrecovered = 0
		corrected = 0
		for pre_barcode in pre_barcode_list:
			corrected_bc = pre_2_post_bc[pre_barcode]
			corrected_bcs_file.write(corrected_bc+'\n')
			post_barcode_list.append(corrected_bc)
			if (pre_barcode != corrected_bc) and corrected_bc:
				corrected +=1
			if not corrected_bc:
				unrecovered += 1
		corrected_bcs_file.close()

		print '\t', total_reads, 'Total reads'
		print '\t', len(pre_barcode_list), "reads are left after TAC-GGG, TSO and polyA filtering" 
		print '\t', len(unique_pre_bcs), 'different barcodes are observed'
		print ''.join(['\t-- ', str(int(round(100*recovered_bc/len(unique_pre_bcs)))), '% of which were recovered'])
		print ''.join(['\t', str(corrected), ' in total, reads were rescued'])		
		print ''.join(['\t', str(int(round(100*(len(pre_barcode_list) - unrecovered)/len(pre_barcode_list)))), '% of reads now have a usable barcode'])
		print "..................................................................."
		os.system("gzip "+f1_base_name+'_processed.fastq')
		print "processed file compressed"

		# read fasta to get gene names
		fasta_return_dicts = get_gene_names_from_fasta(fasta_path)
		dict_gene_names = fasta_return_dicts[0]
		dict_gene_counter = fasta_return_dicts[1]
		gene_counter = len(dict_gene_counter)


		start_time = time.time()
		#Listing *_processed.fastq files:
		preprocessed_files = os.listdir(dir_path_fastqs)
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
		
		read_sam_list = read_sam(sam_file_name)		
		gene_to_umi_bc_dict =read_sam_list[0]
		dict_quality = read_sam_list[1]
		dict_barcode_counter = read_sam_list[2]
		dict_genes_barcode = read_sam_list[3]
		dis_redund = read_sam_list[4]
		read_has_no_bc = read_sam_list[5]
		sam_star = read_sam_list[6]
		bowtie_score = read_sam_list[7]
		dict_quality_scores = read_sam_list[8]

		reads_counted = 0
		bc_2_column = {centroid:i+1 for (i, centroid) in enumerate(centroids)}

		matrix = [[0 for x in range(len(centroids)+1)] for x in range(gene_counter+1)]
		
		# write barcode header for each column
		for centroid in centroids:
			col_num = bc_2_column[centroid]
			matrix[0][col_num] = centroid

		# write the gene name header for each row
		for key_gene in dict_gene_names:
			row_num = dict_gene_counter[dict_gene_names[key_gene]]
			matrix[row_num][0] = dict_gene_names[key_gene]
		
		# write read counts for each gene at each barcode
		for key_gene in dict_gene_names:
			row_num = dict_gene_counter[dict_gene_names[key_gene]]
			for key_barcode in dict_genes_barcode[key_gene]:
				col_num = bc_2_column[key_barcode]
				matrix[row_num][col_num]+=dict_genes_barcode[key_gene][key_barcode]
				reads_counted += dict_genes_barcode[key_gene][key_barcode]
		
		print "Genes-cells matrix created........................................."
		#name=os.path.basename(os.path.normpath(common_path))
		name = f1.split('/')[-1].split('.')[0]
		matrix_file = open(common_path+name+'_matrix.txt', 'w+')
		for item in matrix:
				matrix_file.write('\t'.join([str(i) for i in item])+'\n')
		matrix_file.close()

#### 	call get_html_report
		generate_html_report(common_path, name, bowtie_score, reads_counted, total_reads, preprocessing_saved_reads, 
						 dismissed_reads, dis_tso, dis_no_tac, dis_redund, dict_quality, dict_quality_scores)

		# put read count matrices into list, to save in one matlab file
		matrix_list.append(matrix)
		matrix_name_list.append(name)

		write_to_mat_file(matrix_list, common_path+name)
print "\n"
print "**********************************"
print "**---------Pipeline end---------**"
print "**********************************"
print "\n"