#######################################
## drop_titan_params.py allows you to change
## the parameters for the drop seq exp
## eriment pipeline.
#######################################
## Developed by:
## Paul Rivaud
## paulrivaud.info@gmail.com
## Summer 2015
#######################################


#######################################
## 						      		 ##
## Variables the user might modify:  ##
##				         			 ##
#######################################

####Path to .fastq files
dir_path_fastqs = '../../../RNAseq_files/DS2bis/fastqs/'
#dir_path_fastqs = '../../data/DS2/fastqs/'
#dir_path_fastqs = '../../data/drop_test/fastqs/'
#dir_path_fastqs = '../../data/fastqs_titan/fastqs/'
#dir_path_fastqs = '../../data/fastq_micro/fastqs/'

####Path to .sam files
dir_path_alignment = '../../../RNAseq_files/DS2bis/alignment/'
#dir_path_alignment = '../data/DS2/alignment/'
#dir_path_alignment = '../data/drop_test/alignment/'
#dir_path_alignment = '../data/fastqs_titan/alignment/'
#dir_path_alignment = '../data/fastq_micro/alignment/'

####Path to the reference genome files
reference_genome = '../../reference_genomes/mm9/Transcriptome/transcriptome'
#reference_genome = '../../reference_genomes/Homo_sapiens/UCSC/hg19/cds'
#reference_genome = '../../reference_genomes/Human_Mouse/transcriptome'

#### Path to the species fasta file
fasta_path = '../../reference_genomes/mm9/Transcriptome/transcripts.fa'


#######################################
## 									 ##
## Variables that may not be modifed ##
##(Unless you know what you're doing)##
##									 ##
#######################################

#######################
## Preprocessing var ##
#######################

####TAC length:
tac_length = 3

####umi length:
umi_length = 5

####barcode length
barcode_length = 12

####Tso sequence
tso = 'AAGCAGTGGTATCAACGCAGAGTAC'

####Occurence threshold:
occ_threshold = 5000

#######################
##    Bowtie2 var    ##
#######################
bowtie2_dir=''
#Server processors
processors = 8
#Less sensitive 
bowtie2_options = ['-D', '10', '-R', '1', '-N', '0', '-L', '20','-i', 'S,1,0.50', '-p', str(processors), '--reorder', '--no-hd', '--met', '20']
#very sensitive
#bowtie2_options = ['-D', '20', '-R', '3', '-N', '0', '-L', '20','-i', 'S,1,0.50', '-p', str(processors), '--reorder', '--no-hd', '--met', '20']
