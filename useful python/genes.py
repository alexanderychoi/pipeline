sum_path = os.path.commonprefix([dir_path_fastqs, dir_path_alignment])
if os.path.isfile(sum_path+'genes.txt') and modified == False:
	print "genes.txt already exists... 100%"
else:
	if os.path.isfile(sum_path+'genes.txt'):
		os.remove(sum_path+'genes.txt')
	#Listing .sam files:
	sam_files = os.listdir(dir_path_alignment)
	sam_files.sort()
	nb_sam = len(sam_files)
	curr_sam = 1
	file_genes = open(sum_path+'genes.txt', 'w+')
	print "Creating gene dictionary..."
	gene_dictionary = defaultdict(int)
	print "Checking SAM files..."
	for sam_file in sam_files:
		current_sam = open(dir_path_alignment+sam_file,'r')
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