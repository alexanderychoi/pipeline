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

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

def round_figures(x, n):
	return round(x, int(n - math.ceil(math.log10(abs(x)))))

def to_coord(string,nucleotide_integer,bc_occ):
	coord = "".join(nucleotide_integer.get(base, base) for base in string)
	ints=[]
	ints.extend(coord)
	ints=map(int,ints)
	mult=[bc_occ]*50
	ints=map(operator.mul, ints,mult)	
	return ints

def compare(str1,str2):
	diff = 0
	for i in range(len(str1)):
		if str1[i]!=str2[i]:
			diff+=1
	return diff

def bc_dist(s1, s2):
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

def measure_dist(bc, centroids_list):
	dist_to_centroids=[0 for i in range(len(centroids_list))]
	for i in range(len(centroids_list)):
		dist_to_centroids[i]=bc_dist(bc,centroids_list[i])
	return dist_to_centroids

def cluster_points(bc_list, centroids_list):
	clusters  = [[] for i in range(len(centroids_list))]
	bc_2_centroid=defaultdict(str)
	for bc in bc_list:
		dist_to_centroids = measure_dist(bc, centroids_list)
		mini=min(dist_to_centroids)
		ind=[i for i, j in enumerate(dist_to_centroids) if j == mini]
		clusters[ind[0]].append(bc)
		bc_2_centroid[bc]=centroids_list[ind[0]]
	return [clusters,bc_2_centroid]

def update_centroids(centroids_list,clusters):
	new_centroid_list=centroids_list
	for i in range(len(clusters)):
		mean_coords=[0]*50
		number_of_bc=0
		for bc in clusters[0][i]:
			bc_occ=barcode_count_dict[bc]
			number_of_bc+=bc_occ
			coords=to_coord(bc,nucleotide_integer,bc_occ)
			mean_coords=map(operator.add, mean_coords,coords)
		if number_of_bc == 0:
			number_of_bc_list=[1]*50
		else:
			number_of_bc_list=[float(number_of_bc)]*50
		mean_coords=map(operator.div,mean_coords,number_of_bc_list)
		new_centroid=''
		low_ind=0
		high_ind=5
		while high_ind<=50:
			sub=[]
			for item in mean_coords[low_ind:high_ind]:
				sub.append(item)
			m=max(sub)
			ind=[i for i, j in enumerate(sub) if j == m]
			if ind[0]==0:
				new_nucleotide='A'
			elif ind[0]==1:
				new_nucleotide='G'
			elif ind[0]==2:
				new_nucleotide='C'
			elif ind[0]==3:
				new_nucleotide='T'
			elif ind[0]==4:
				new_nucleotide='N'
			new_centroid+=new_nucleotide
			low_ind+=5
			high_ind+=5
		new_centroid_list[i]=new_centroid
	return new_centroid_list

def cluster_pop_to_string(my_clusters):
	# print length of clusters
	for c in range(len(my_clusters)):
		print "c"+str(c)+": "+str(len(my_clusters[c]))

def population_to_list(iteration):
	pop_list = [[] for i in range(len(iteration))]
	for i in range(len(iteration)):
		pop_list[i] = len(iteration[i])
	return pop_list

def kmeans(bc_list,centroids_list,previous_population):
	iteration = cluster_points(bc_list, centroids_list)
	new_population = population_to_list(iteration)
	print new_population
	if previous_population != new_population:
		updated_centroids = update_centroids(centroids_list,iteration)
		print updated_centroids
		kmeans(bc_list,updated_centroids,new_population)
	else:
		print "End of K-MEANS"
		print "Returning object [centroids (list), clusters (list of lists)]"
		return [centroids_list,iteration]

def kmeans_one_iter(bc_list,centroids_list):
	iteration = cluster_points(bc_list, centroids_list)
	updated_centroids = update_centroids(centroids_list,iteration)
	return [iteration, updated_centroids]


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
barcode_list = []
barcode_count_dict = defaultdict(int)

for f1, f2 in grouped(fastq_files, 2):
	percent = int((curr_fastq/nb_fastqs)*100)
	curr_fastq+=2
	f1_base_name = f1.split(os.extsep)[0]

	if os.path.isfile(dir_path_fastqs+f1.split(os.extsep)[0]+'_processed.fastq.gz'):
		print f1
		print f2

		pre_bc_file = open(f1_base_name + '_pre_barcode.txt', 'r', 1)
		barcode_list = pre_bc_file.read().splitlines()
		pre_bc_file.close()

		umi_file = open(f1_base_name + '_umi.txt', 'r', 1)
		umi_list = umi_file.read().splitlines()
		umi_file.close()


		print "\tFiles already preprocessed. Loading in parameters ...",percent,"%"
	else:
		print "Begin..."
		print "Opening fastq files..."
		file1_fastq = gzip.open(f1,'rb')
		file2_fastq = gzip.open(f2,'rb')
		print f1
		print f2
		print "\tReading files..."
		print "\tWriting files..."
		file_processed = open(f1_base_name+'_processed.fastq', 'w+', 1)
		file_umi = open(f1_base_name+'_umi.txt', 'w+',1)
		file_pre_barcode = open(f1_base_name+'_pre_barcode.txt', 'w+', 1)
		#Stats about the trimming process
		total_reads = 0
		saved_reads = 0
		dismissed_reads = 0
		dis_tso = 0
		dis_no_tac = 0
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
			if (str_search in f1_line2[:3]) and f1_line2[:6] != (str_search + 'GGG'):
				if tso not in f2_line2:
					barcode = f1_line2[tac_length:tac_length+barcode_length]
					umi = f1_line2[tac_length+barcode_length:tac_length+barcode_length+umi_length]
					file_umi.write(umi+'\n')
					umi_list.append(umi)
					barcode_list.append(barcode)					
					barcode_count_dict[barcode] += 1
					file_pre_barcode.write(barcode+'\n')

					# write read (4 lines) to processed file 
					file_processed.write(f1_line1 + f2_line2 + f1_line3 + f2_line4)

					saved_reads+=1
				else:
					dis_tso+=1
					dismissed_reads+=1
			else:
				dis_no_tac+=1
				dismissed_reads+=1

		file_pre_barcode.close()
		file1_fastq.close()
		file2_fastq.close()
		file_processed.close()


		#################################################################
		#################################################################
		#################################################################
		#################################################################
		##### KMEANS ####################################################
		#################################################################
		#################################################################
		#################################################################
		#################################################################

		#pre_bc_file = open(f1_base_name+'_pre_barcode.txt', 'r', 1)
		nucleotide_integer = {'A': '10000', 'T': '00010', 'C': '00100', 'G': '01000', 'N': '00001'}
		

		# sort barcodes on occurances to get top 45
		dd = sorted(barcode_count_dict.items(), key=operator.itemgetter(1))
		short_list={}
		for i in range(len(dd)):
			if i>=(len(dd)-45):
				short_list[dd[i][0]] = dd[i][1]

		l=sorted(short_list,key=short_list.get)
		l=l[::-1]
		sim=[]
		for i in range(len(l)):
			for j in range(len(l)):
				if j>i:
					diff=bc_dist(l[i],l[j])
					if diff<=2:
						sub=[l[i],l[j]]
						i_in_sim=0
						i_value_for_x=0
						j_in_sim=0
						j_value_for_x=0
						for x in range(len(sim)):
							if sub[0] in sim[x]:
								i_in_sim=1
								i_value_for_x=x
							if sub[1] in sim[x]:
								j_in_sim=1
								j_value_for_x=x
						if i_in_sim==1 and j_in_sim==1:
							sim[i_value_for_x]=sim[i_value_for_x]+sim[j_value_for_x]
							sim.pop(j_value_for_x)
						elif i_in_sim==1 and j_in_sim==0:
							sim[i_value_for_x].append(sub[1])
						elif i_in_sim==0 and j_in_sim==1:
							sim[j_value_for_x].append(sub[0])
						else:
							sim.append(sub)
					else:
						i_in_sim=0
						j_in_sim=0
						for y in range(len(sim)):
							if l[i] in sim[y]:
								i_in_sim=1
							if l[j] in sim[y]:
								j_in_sim=1
						if i_in_sim==0:
							sim.append([l[i]])
						if j_in_sim==0:
							sim.append([l[j]])
		centroids={}
		value=0
		for k in range(len(sim)):
			value=0
			for m in range(len(sim[k])):
				value+=short_list[sim[k][m]]
			centroids[sim[k][0]]=value

		centroids = dict(sorted(centroids.iteritems(), key=operator.itemgetter(1), reverse=True)[:cell_num])
		centroids_list=[]
		for key in centroids:
			centroids_list.append(key)
			
		print len(barcode_list)
		#begin_population=[0]*15
		#kmeans(dict_bc.keys(),centroids_list,begin_population)
		res=kmeans_one_iter(barcode_list,centroids_list)
		saved=0
		for i in range(len(res[0][0])):
			saved+=len(res[0][0][i])
		print saved
		#read pre barcode file
		pre_bc=open(f1_base_name+'_pre_barcode.txt', 'r', 1)
		clustered_bc=open(f1_base_name+'_barcode.txt', 'w+', 1)
		corrected_barcode_list = []
		for barcode in barcode_list:
			if barcode in res[0][1].keys():
				clustered_bc.write(res[0][1][barcode]+'\n')
				corrected_barcode_list.append(res[0][1][barcode])
			else:
				corrected_barcode_list.append(barcode)
				clustered_bc.write(barcode+'\n')
		clustered_bc.close()
		pre_bc.close()




		print "\tTotal reads: ",total_reads
		print "...................................................................",percent,"%"
		os.system("gzip "+f1_base_name+'_processed.fastq')
		print "processed file compressed"

		############################### ALIGNMENT ###############################
		start_time = time.time()
		#Listing *_processed.fastq files:
		preprocessed_files = os.listdir(dir_path_fastqs)
		#Cleaning list of files in the directory (Only *processed.fastq files):
		files_processed = [f for f in preprocessed_files if 'processed.fastq.gz' in f]
		files_processed.sort()
		bowtie_opt = ' '.join(bowtie2_options)
		print "\n"
		print "**********************************"
		print "**      Starting alignment      **"
		print "**********************************"
		print "\n"

		nb_processed = len(files_processed)
		curr_processed = 1

		for file_processed in files_processed:
			percent = int((curr_processed/nb_processed)*100)
			curr_processed+=1
			if os.path.isfile(dir_path_alignment+file_processed.split(os.extsep)[0]+'.sam.gz'):
				print "Skipping alignment: "+file_processed.split(os.extsep)[0]+'.sam.gz already exists...',percent,"%"
			else:
				print "Starting alignment with bowtie2 for", file_processed
				os.system("nice "+bowtie2_dir+" bowtie2 "+bowtie_opt+\
					" -x "+reference_genome+" -U "+dir_path_fastqs+file_processed+\
					" -S "+dir_path_alignment+file_processed.split(os.extsep)[0]+'.sam')
				print file_processed+" aligned..."
				print file_processed.split(os.extsep)[0]+'.sam created...'
				print "...................................................................",percent,"%"
				os.system("gzip "+dir_path_alignment+file_processed.split(os.extsep)[0]+'.sam')
				print "Sam file compressed..."
				print "...................................................................",percent,"%"
		total_time = time.time() - start_time
		print "Reads alignment time:"
		print int(total_time/60),"min",int(total_time%60),"sec"

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
			gene_counter = 1
			barcode_counter = 1
			sam_file = gzip.open(dir_path_alignment+file_processed.split(os.extsep)[0]+'.sam.gz','rb')
			list_f = os.listdir(dir_path_fastqs)
			#dict_genes_barcode = defaultdict(dict)
			dict_genes_barcode = defaultdict(lambda : defaultdict(int))
			dict_gene_counter = defaultdict(int)
			dict_gene_names = defaultdict(str)
			dict_barcode_counter = defaultdict(int)
			dict_barcode_occurences = defaultdict(int)

			print "Storing data in dictionaries..."


			###########################################################
			# Get the list of all genes from the genome fasta file
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
			###########################################################

			# collect alignment statistics
			sam_gene=0
			sam_star=0
			bowtie_al=0
			#dict_quality=defaultdict(dict)
			dict_quality = defaultdict(lambda : defaultdict(int))
			dict_quality_scores=defaultdict(int)
			sam_line_ind = 0
			dis_redund = 0

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
						print gene, sam_line_ind
						bowtie_al+=1
						umi = umi_list[sam_line_ind]
						barcode = corrected_barcode_list[sam_line_ind]
						umi_bc_pair = (umi, barcode)
						if umi_bc_pair not in gene_to_umi_bc_dict[gene]:
							print 'about to check quality score'
							AS_score = int(columns[11][5:])
							AS_score_str = str(AS_score)
							dict_quality_scores[AS_score_str]+=1
							if AS_score>=-3:
								print gene, ' passed quality check'
								sam_gene+=1
								gene_to_umi_bc_dict[gene].add(umi_bc_pair)
								dict_quality[barcode]['high']+=1
								dict_barcode_counter[barcode] +=1
								barcode_counter+=1
								dict_genes_barcode[gene][barcode] += 1
							else:
								dict_quality[barcode]['low']+=1
								sam_star += 1
						else:
							dis_redund += 1
					else:
						sam_star += 1
				sam_line_ind += 1

			# compute read count and alignment stats	
			alignment_score=sam_gene/(sam_gene+sam_star)
			alignment_score=round(alignment_score*100)
			bowtie_score=bowtie_al/(bowtie_al+sam_star)
			bowtie_score=round(bowtie_score*100)
			sam_file.close()
			print "Data stored in dictionaries........................................",percent,"%"
			print "Creating genes-cells matrix...\n"
			print cell_num, "cells"			
			print gene_counter-1, "genes"
			print barcode_counter-1, "counted reads"

			bc_2_column = {key:(value+1) for (value, key) in enumerate(centroids.keys())}

			matrix = [[0 for x in range(cell_num+1)] for x in range(gene_counter)]
			
			# write barcode header for each column
			for key_barcode in dict_barcode_counter:
				col_num = bc_2_column[key_barcode]
				matrix[0][col_num] = key_barcode

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
			
			print "Genes-cells matrix created.........................................",percent,"%"
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
			html_file.write(str(alignment_score))
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
			html_file.write(str(saved_reads))
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