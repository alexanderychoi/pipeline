# -*-coding:Utf-8 -*

import shutil
from collections import defaultdict

matrix1_path="./sample2_truncBarcode_matrix.txt"
matrix1stats_path="./sample2_truncBarcode_matrix_stats.txt"

human_path="../../reference_genomes/Homo_sapiens/UCSC/hg19/cds.fa"
mouse_path="../../reference_genomes/mm9/Transcriptome/transcripts.fa"

human_file = open(human_path,'r')
set_human_genes = defaultdict(str)

while True:
	line = human_file.readline()
	if not line:
		break
	else:
		if line[:1] == '>':
			gene_NM = line.split(' ')[0]
			gene_NM=gene_NM[1:]
			gene_NM=gene_NM.replace('\n','')
			gene_NM=gene_NM.upper()

			gene_real_name = line.split(' ')[1]
			gene_real_name=gene_real_name[5:]
			gene_real_name=gene_real_name.replace("\n","")
			gene_real_name=gene_real_name.upper()
			set_human_genes[gene_NM.replace(" ","")]=gene_real_name.replace(" ","")

mouse_file = open(mouse_path,'r')
set_mouse_genes = defaultdict(str)
while True:
	line = mouse_file.readline()
	if not line:
		break
	else:
		if line[:1] == '>':
			gene_NM = line.split(' ')[0]
			gene_NM=gene_NM[1:]
			gene_NM=gene_NM.replace('\n','')
			gene_NM=gene_NM.upper()

			gene_real_name = line.split(' ')[1]
			gene_real_name=gene_real_name[5:]
			gene_real_name=gene_real_name.replace("\n","")
			gene_real_name=gene_real_name.upper()
			set_mouse_genes[gene_NM.replace(" ","")]=gene_real_name.replace(" ","")

mouse_file.close()
human_file.close()
###############################################

matrix1 = open(matrix1_path,'r')
matrix1_stats = open(matrix1stats_path,'w+')

line=matrix1.readline()
columns = line.split("\t",1)
matrix1_stats.write(columns[0]+"\tM/H\t"+columns[1])

while True:
	line=matrix1.readline()
	if not line:
		break
	else:
		columns = line.split("\t",1)
		gene = columns[0]
		gene=gene.replace(" ","")
		if gene in set_mouse_genes:
			matrix1_stats.write(set_mouse_genes[gene]+"\t0\t"+columns[1])
		elif gene in set_human_genes:
			matrix1_stats.write(set_human_genes[gene]+"\t1\t"+columns[1])

matrix1.close()
matrix1_stats.close()
matrix1_stats = open(matrix1stats_path,'r+')

line = matrix1_stats.readline()
columns = line.split('\t')
list_lenght = len(columns)-2 
human_reads = [0]*(list_lenght)
total_reads = [0]*(list_lenght)
specificity = [0]*(list_lenght)
while True:
	line = matrix1_stats.readline()
	if not line:
		break
	else:
		columns = line.split('\t')
		for i in range(2,list_lenght+2):
			total_reads[i-2] += int(columns[i])
			human_reads[i-2] += int(columns[i])*int(columns[1])

for i in range(0,list_lenght):
	human_specificity = round(human_reads[i]*100/total_reads[i],1)
	mouse_specificity = 100-human_specificity
	if human_specificity>mouse_specificity:
		specificity[i]=str(human_specificity)+"-human"
	else:
		specificity[i]=str(mouse_specificity)+"-mouse"

matrix1_stats.write("0\tHuman Reads\t"+'\t'.join([str(i) for i in human_reads])+'\n')
mouse_reads = [int(j)-int(i) for i,j in zip(human_reads,total_reads)]
matrix1_stats.write("0\tMouse Reads\t"+'\t'.join([str(i) for i in mouse_reads])+'\n')
matrix1_stats.write("0\tTotal Reads\t"+'\t'.join([str(i) for i in total_reads])+'\n')
matrix1_stats.write("0\tSpecificity\t"+'\t'.join([str(i) for i in specificity])+'\n')


matrix1_stats.close()


