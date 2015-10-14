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
#dir_path_fastqs = '/home/graham/dropseq_pipe/drop_test_files/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/pipeline/targeted_test_files/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/targeted_test_files/fastqs/'
dir_path_fastqs = '/home/graham/dropseq_pipe/DS8/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/miseq_test/fastqs/'
#dir_path_fastqs = '/home/david/RNAseq_files/DS4-H2/fastqs/'
#dir_path_fastqs = '/home/david/RNAseq_files/DS8/fastqs/'

####Path to .sam files
#dir_path_alignment = '/home/graham/dropseq_pipe/drop_test_files/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/pipeline/targeted_test_files/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/targeted_test_files/alignment/'
dir_path_alignment = '/home/graham/dropseq_pipe/DS8/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/miseq_test/alignment/'
#dir_path_alignment = '/home/david/RNAseq_files/DS4-H2/alignments/'
#dir_path_alignment = '/home/david/RNAseq_files/DS8/alignments/'

####Path to the reference genome files
#reference_genome = '~/reference_genomes/mm9/Transcriptome/transcriptome'
reference_genome = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcriptome'

#reference_genome = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcriptome'

#### Path to the species fasta file
fasta_path = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcripts.fa'
#fasta_path = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcripts.fa'


#######################################
## 									 ##
## Variables that may not be modifed ##
##(Unless you know what you're doing)##
##									 ##
#######################################

#######################
## Preprocessing var ##
#######################

#### String to search at beginning of R1:
str_search='TAC'
#str_search=''

#### umi length:
umi_length = 9

#### barcode length
barcode_length = 10

####Tso sequence
tso = 'AAGCAGTGGTATCAACGCAGAGTAC'

#### Number of cells expected (same as # of clusters)
cell_num = 5


#### threshold for filtering consecutive A's
A_num = 8

#### Number of bases required for a read, throws out anything shorter
min_read_len = 30

#### Occurence threshold:
#occ_threshold = 5000

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



############################################################################################
#######################         !! DONT TOUCH !!        ####################################
############################################################################################

# compute length of tac for use as constant
tac_length = len(str_search)

# what length of seq is considered a TSO match?
tso_word_len = 8




############################################################################################
############################################################################################
############################################################################################
