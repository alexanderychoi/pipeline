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
import re
import math
import os
import os.path
import time
import gzip
import operator
import glob
import sklearn
import numpy as np
from matplotlib import pyplot as plt
from scipy import cluster, spatial
import scipy

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def compare(str1,str2):
	diff = 0
	for i in range(len(str1)):
		if str1[i]!=str2[i]:
			diff+=1
	return diff

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
				gene_symbol=gene_symbol.strip()
				gene_symbol=gene_symbol.upper()
				if gene_symbol not in dict_gene_counter:
					dict_gene_counter[gene_symbol] = gene_counter
					gene_counter+=1
				gene_nm=line.split(' ')[0]
				gene_nm=gene_nm[1:]
				gene_nm=gene_nm.strip()
				gene_nm=gene_nm.upper()
				dict_gene_names[gene_nm] = gene_symbol
	fasta_file.close()
	return [dict_gene_names, dict_gene_counter]

def preprocess_fastq_file_pair(f1, f2, str_search):
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
	while True:
		f1_line1 = file1_fastq.readline()
		f1_line2 = file1_fastq.readline()
		f1_line3 = file1_fastq.readline()
		f1_line4 = file1_fastq.readline()
		f2_line1 = file2_fastq.readline()
		f2_line2 = file2_fastq.readline()
		f2_line3 = file2_fastq.readline()
		f2_line4 = file2_fastq.readline()
		if not f1_line1:
			break
		else:
			total_reads+=1
		if (str_search in f1_line2[:tac_length]) and f1_line2[:(tac_length+3)] != (str_search + 'GGG'):
			if tso not in f2_line2:
				barcode = f1_line2[tac_length:tac_length+barcode_length]
				umi = f1_line2[tac_length+barcode_length:tac_length+barcode_length+umi_length]
				file_umi.write(umi+'\n')
				umi_list.append(umi)
				pre_barcode_list.append(barcode)					
				barcode_count_dict[barcode] += 1
				file_pre_barcode.write(barcode+'\n')

				# write read (4 lines) to processed file 
				file_processed.write(f1_line1 + f2_line2 + f1_line3 + f2_line4)

				preprocessing_saved_reads+=1
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
		preprocess_results = preprocess_fastq_file_pair(f1, f2, str_search)
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
		top_bc_num = 1500;
		sorted_bcs = sorted(barcode_count_dict.items(), key=operator.itemgetter(1))
		top_bcs={}
		for i in range(max(num_unique_pre_bcs-top_bc_num, 0), num_unique_pre_bcs):
			top_bcs[sorted_bcs[i][0]] = sorted_bcs[i][1]

		bc_dist_mat = compute_bc_dist_mat(top_bcs.keys())
		keep_dist_mat_inds = np.sum(bc_dist_mat<3, axis=0)>(round(top_bc_num*.0035))
		pruned_bc_dist_mat = bc_dist_mat[keep_dist_mat_inds , :]
		pruned_bc_dist_mat = pruned_bc_dist_mat[:, keep_dist_mat_inds]
		pruned_top_bcs = [bc for i,bc in enumerate(top_bcs.keys()) if keep_dist_mat_inds[i]]





		### plot dendrogram and clustered barcodes ###
		# fig = pylab.figure(figsize=(8,8))
		# ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
		Y1 = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.squareform(pruned_bc_dist_mat), method='average')
		# Z1 = scipy.cluster.hierarchy.dendrogram(Y1, orientation='right')
		# ax1.set_xticks([])
		# ax1.set_yticks([])

		# fig = pylab.figure(figsize=(8,8))
		# ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
		# Y2 = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.squareform(pruned_bc_dist_mat), method='average')
		# Z2 = scipy.cluster.hierarchy.dendrogram(Y2)
		# ax2.set_xticks([])
		# ax2.set_yticks([])

		# # Plot distance matrix.
		# axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
		# idx1 = Z1['leaves']
		# idx2 = Z2['leaves']
		# bc_dist_mat = pruned_bc_dist_mat[idx1,:]
		# bc_dist_mat = pruned_bc_dist_mat[:,idx2]
		# im = axmatrix.matshow(pruned_bc_dist_mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
		# axmatrix.set_xticks([])
		# axmatrix.set_yticks([])

		# # Plot colorbar.
		# axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
		# pylab.colorbar(im, cax=axcolor)
		# fig.show()
		#############################################
		clusters = scipy.cluster.hierarchy.fcluster(Y1, 4, criterion='distance')
		num_clust = np.max(clusters)

		# find the bc that has the lowest distance to all other bc's within a cluster
		centroids = []
		for i in range(1,num_clust):
			cluster_dist_mat = pruned_bc_dist_mat[clusters==i,:]
			cluster_dist_mat = cluster_dist_mat[:,clusters==i]
			cluster_bcs = [bc for j,bc in enumerate(pruned_top_bcs) if clusters[j]==i]
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
		for pre_barcode in pre_barcode_list:
			corrected_bc = pre_2_post_bc[pre_barcode]
			corrected_bcs_file.write(corrected_bc+'\n')
			post_barcode_list.append(corrected_bc)
			if not corrected_bc:
				unrecovered += 1
		corrected_bcs_file.close()

		print '\t', len(pre_barcode_list), 'total usable reads observed' 
		print '\t', len(unique_pre_bcs), 'different barcodes'
		print ''.join(['\t', str(int(round(100*recovered_bc/len(unique_pre_bcs)))), '% unique barcodes recovered'])
		print ''.join(['\t', str(int(round(100*(len(pre_barcode_list) - unrecovered)/len(pre_barcode_list)))), '% total barcodes recovered'])
		print '\t', total_reads, 'Total reads'
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
		

		############################### Data gathering ###############################
		print "\n"
		print "**********************************"
		print "**       Retrieving data        **"
		print "**********************************"
		print "\n"

		sum_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])
		if os.path.isfile(sum_path+'matrix.txt'):
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
								#	sam_gene+=1
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
			reads_counted = 0
			sam_file.close()
			print "Data stored in dictionaries........................................"
			print "Creating genes-cells matrix...\n"
			print gene_counter-1, "genes"
			print barcode_counter-1, "counted reads"

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
			name=os.path.basename(os.path.normpath(sum_path))
			matrix_file = open(sum_path+name+'_matrix.txt', 'w+')
			for item in matrix:
					matrix_file.write('\t'.join([str(i) for i in item])+'\n')
			matrix_file.close()


			html_file = open(sum_path+name+'_report.html', 'w+')
			html_file.write('''<!DOCTYPE html>
			<html>
				<head>
					<link rel="stylesheet" type="text/css" href="https://bootswatch.com/cerulean/bootstrap.min.css">
				</head>
				<body>
					<nav class="navbar navbar-default">
					  <div class="container-fluid">
					    <div class="navbar-header">
					      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#bs-example-navbar-collapse-1">
					        <span class="sr-only">Toggle navigation</span>
					        <span class="icon-bar"></span>
					        <span class="icon-bar"></span>
					        <span class="icon-bar"></span>
					      </button>
					      <a class="navbar-brand" onclick="return false">HTML Report for ''')
			html_file.write(name)
			html_file.write('''</a>
					    </div>
					    <div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
					      	<ul class="nav navbar-nav navbar-right">
				        		<li><a onclick="return false">Thomson Lab</a></li>
				      		</ul>
					    </div>
					  </div>
					</nav>
					<br>
					<h1 class="text-primary">Preprocessing results</h1>
					<br>
					<div id="preprocessing_bar" style="height: 300px; width: 50%;"></div>
					<br>
					<h1 class="text-primary">Dismissed reads information</h1>
					<br>
					<div id="dismissed_info" style="height: 300px; width: 50%;"></div>
					<br>
					<h1 class="text-primary">Bowtie2 overall alignment rate</h1>
					<br>
					<h1>''')
			html_file.write(str(bowtie_score))
			html_file.write('''%</h1>
					<br>
					<h1 class="text-primary">Reads used for gene expression quantification</h1>
					<br>
					<h1>''')
			html_file.write(str(100*round(reads_counted/total_reads,3)))
			html_file.write('''%</h1>
					<br>
					<h1 class="text-primary">Reads per cell distribution</h1>
					<br>
					<div id="readsDistribution" style="height: 300px; width: 50%;"></div>
					<br>
					<h1 class="text-primary">Quality scores histogram</h1>
					<br>
					<div id="qualityScore" style="height: 300px; width: 50%;"></div>
					


					</body>
				</html>


				<script type="text/javascript" src="http://canvasjs.com/assets/script/canvasjs.min.js"></script>
				<script type="text/javascript">
					function preprocessingPlot() {
						var chart = new CanvasJS.Chart("preprocessing_bar", {
									theme: "theme2",//theme1
									title:{
										text: ""
									},
									animationEnabled: true,   // change to true
									data: [              
									{
										// Change type to "bar", "splineArea", "area", "spline", "pie",etc.
										type: "column",
										dataPoints: [
											{ label: "Total reads",  y: ''')
			html_file.write(str(total_reads))
			html_file.write('''  },
								{ label: "Saved reads", y: ''')
			html_file.write(str(preprocessing_saved_reads))
			html_file.write('''  },
								{ label: "Dismissed reads", y: ''')
			html_file.write(str(dismissed_reads))
			html_file.write('''  }
										]
									}
									]
								});
								chart.render();
					}

					function dismissedInfo(){
						var chart1 = new CanvasJS.Chart("dismissed_info",
									{
										title:{
											text: ""
										},
							                        animationEnabled: true,
										theme: "theme2",
										data: [
										{        
											type: "doughnut",
											indexLabelFontFamily: "Garamond",       
											indexLabelFontSize: 20,
											startAngle:0,
											indexLabelFontColor: "dimgrey",       
											indexLabelLineColor: "darkgrey", 
							

											dataPoints: [
											{  y: ''')
			html_file.write(str(dis_tso))
			html_file.write(''', label: "Contains TSO" },
								{  y: ''')
			html_file.write(str(dis_no_tac))
			html_file.write(''', label: "No TAC" },
								{  y: ''')
			html_file.write(str(dis_redund))
			html_file.write(''', label: "Redundant" }

											]
										}
										]
									});
									chart1.render();
					}

					function distribution(){
						var chart2 = new CanvasJS.Chart("readsDistribution",
						    {
						      title:{
						      text: ""   
						      },
						      axisY:{
						        title:"Number of reads"   
						      },
						      animationEnabled: true,
						      data: [
						      {        
						        type: "stackedColumn",
						        toolTipContent: "{label}<br/><span style='\\"'color: {color};'\\"'><strong>{name}</strong></span>: {y} reads",
						        name: "Low quality",
						        showInLegend: "true",
						        dataPoints: [''')
			for key in sorted(dict_quality.keys()):
				html_file.write('{  y: ')
				html_file.write(str(dict_quality[key]['low']))
				html_file.write(', label:"')
				html_file.write(key)
				html_file.write('"},\n')
			html_file.write(''']

						      },  {        
						        type: "stackedColumn",
						        toolTipContent: "{label}<br/><span style='\\"'color: {color};'\\"'><strong>{name}</strong></span>: {y} reads",
						        name: "Good quality",
						        showInLegend: "true",
						        dataPoints: [''')
			for key in sorted(dict_quality.keys()):
				html_file.write('{  y: ')
				html_file.write(str(dict_quality[key]['high']))
				html_file.write(', label:"')
				html_file.write(key)
				html_file.write('"},\n')
			html_file.write(''']
						      }            
						      ]
						      ,
						      legend:{
						        cursor:"pointer",
						        itemclick: function(e) {
						          if (typeof (e.dataSeries.visible) ===  "undefined" || e.dataSeries.visible) {
							          e.dataSeries.visible = false;
						          }
						          else
						          {
						            e.dataSeries.visible = true;
						          }
						          chart2.render();
						        }
						      }
						    });

						    chart2.render();
					}

					function quality(){
						var chart3 = new CanvasJS.Chart("qualityScore",
							{
								animationEnabled: true,
								title:{
									text: ""
								},
								data: [
								{
									type: "column", //change type to bar, line, area, pie, etc
									dataPoints: [''')
			for key in dict_quality_scores.keys():
				html_file.write('{  x: ')
				html_file.write(key)
				html_file.write(', y: Math.log10(')
				html_file.write(str(dict_quality_scores[key]))
				html_file.write(')},\n')
			html_file.write(''']
								}
								]
							});

							chart3.render();
					}
					function loadAll() {
						preprocessingPlot();
						dismissedInfo();
						distribution();
						quality();
					}

					window.onload = loadAll;
				</script>''')
			html_file.close()

		print "HTML report completed"

		print "\n"
		print "**********************************"
		print "**---------Pipeline end---------**"
		print "**********************************"
		print "\n"