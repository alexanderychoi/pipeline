######################################
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
#dir_path_fastqs = '/home/graham/dropseq_pipe/pipeline/targeted_dropseq_test/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/pipeline/targeted_dropseq_test/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/DS9/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/miseq_test/fastqs/'
#dir_path_fastqs = '/home/david/RNAseq_files/DS4-H2/fastqs/'
#dir_path_fastqs = '/home/david/RNAseq_files/DS8/fastqs/'
#dir_path_fastqs = '/home/david/TITAN_pipeline/pipeline_code/targeted_dropseq_test/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/pipeline/sisi_dropseq_test/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/DS4_cc/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/NativeDiff/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/DS22/very_sensitive/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/DS6-H/fastqs/'
#dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/DS5/fastqs/'
dir_path_fastqs = '/home/iamcam/Documents/RNAseqdata/DS35/fastqs/'

####Path to .sam files
#dir_path_alignment = '/home/graham/dropseq_pipe/drop_test_files/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/pipeline/targeted_dropseq_test/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/pipeline/targeted_dropseq_test/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/DS9/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/miseq_test/alignment/'
#dir_path_alignment = '/home/david/RNAseq_files/DS4-H2/alignments/'
#dir_path_alignment = '/home/david/RNAseq_files/DS8/alignments/'
#dir_path_alignment = '/home/david/TITAN_pipeline/pipeline_code/targeted_dropseq_test/alignment/'
#dir_path_alignment = '/home/iamcam/Documents/pipeline/sisi_dropseq_test/alignment/'
#dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/DS4_cc/alignments/'
#dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/NativeDiff/alignments/'
#dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/DS22/very_sensitive/alignments/'
#dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/DS6-H/alignments/'
#dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/DS5/alignments/'
dir_path_alignment = '/home/iamcam/Documents/RNAseqdata/DS35/alignments/'

####Path to the reference genome files
#reference_genome = '~/reference_genomes/mm9/Transcriptome/transcriptome'
#reference_genome = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcriptome'
#reference_genome = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcriptome'
#reference_genome = '/home/iamcam/Documents/reference_genomes/human_mouse/human_mouse'
reference_genome = '/home/iamcam/Documents/reference_genomes/mm9/Transcriptome/mouse_cds'
#reference_genome = '/home/iamcam/Documents/reference_genomes/barcodes/barcodes'
#reference_genome = '/home/iamcam/Documents/reference_genomes/Homo_sapiens/UCSC/hg19/human_cds'


#### Path to the species fasta file
#fasta_path = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcripts.fa'
#fasta_path = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcripts.fa'
fasta_path = '/home/iamcam/Documents/reference_genomes/mm9/Transcriptome/mouse_transcripts.fa'
#fasta_path = '/home/iamcam/Documents/reference_genomes/barcodes/barcodes.fa'
#fasta_path = '/home/iamcam/Documents/reference_genomes/Homo_sapiens/UCSC/hg19/human_cds.fa'

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
#str_search='TAC'
str_search=''

#### umi length:
umi_length = 8

#### barcode length
barcode_length = 12

####Tso sequence
tso = 'AAGCAGTGGTATCAACGCAGAGTAC'

#### Number of cells expected (same as # of clusters)
cell_num = 50

#### threshold for filtering consecutive A's
A_num = 8

#### Number of bases required for a read, throws out anything shorter
min_read_len = 30

#number of frame shift bases allowed
fs =  6 #number of frame shift bases allowed

#### Occurence threshold:
#occ_threshold = 5000

#######################
##    Bowtie2 var    ##
#######################

bowtie2_dir=''
#Server processors
processors = 8
#Less sensitive 
#bowtie2_options = ['-D', '10', '-R', '1', '-N', '0', '-L', '20','-i', 'S,1,0.50', '-p', str(processors), '--reorder', '--no-hd', '--met', '20']
#very sensitive
bowtie2_options = ['-D', '20', '-R', '3', '-N', '0', '-L', '20','-i', 'S,1,0.50', '-p', str(processors), '--reorder', '--no-hd', '--met', '20']

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
