#######################################
## plate_params.py allows you to change
## the parameters for the plate experim
## ent pipeline.
#######################################
## Developed by:
## Paul Rivaud
## paulrivaud.info@gmail.com
## Summer 2015
#######################################


#######################################
## 									 ##
## Variables the user might modify:  ##
##									 ##
#######################################

####Path to .fastq files
#dir_path_fastqs = '../../data/fastqs_lucas/fastqs/'
#dir_path_fastqs = '/home/graham/dropseq_pipe/pipeline/bulk_test/fastqs/'
dir_path_fastqs = '/home/david/RNAseq_files/NeuralDiff_bulk/fastqs/'

####Path to .sam files
#dir_path_alignment = '../../data/fastqs_lucas/alignment/'
#dir_path_alignment = '/home/graham/dropseq_pipe/pipeline/bulk_test/alignment/'
#dir_path_alignment = '../data/test/alignment/'
dir_path_alignment = '/home/david/RNAseq_files/NeuralDiff_bulk/alignments/'

####Path to the reference genome files
#reference_genome = '../reference_genomes/mm9/Transcriptome/transcriptome'
#reference_genome = '../../reference_genomes/Rattus_norvegicus/UCSC/rn5/cds'
#reference_genome = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcriptome'
reference_genome = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcriptome'

#### Path to the species fasta file
#fasta_path = '../../reference_genomes/Rattus_norvegicus/UCSC/rn5/cds.fa'
#fasta_path = '/home/graham/dropseq_pipe/reference_genomes/mm9/Transcriptome/transcripts.fa'
fasta_path = '/home/david/TITAN_pipeline/reference_genomes/mm9/Transcriptome/transcripts.fa'

#######################################
## 									 ##
## Variables that may not be modifed ##
##(Unless you know what you're doing)##
##									 ##
#######################################

#######################
## Preprocessing var ##
#######################

####umi length:
umi_length = 0

####TAC length:
tac_length = 0

####Number of G's to delete
g_length = 0

####Tso sequence
tso = 'AAGCAGTGGTATCAACGCAGAGTAC'

####Number of Ts or As in a row beyond which the polyT or polyA sequence is trimmed
num_ta = 10
error_string_t = 'T'*num_ta
error_string_a = 'A'*num_ta

####Minimum length for RNA seq (if < minimum_length,  the read is dismisse)
minimum_length = 30


#######################
##    Bowtie2 var    ##
#######################
bowtie2_dir=''
#Server processors
processors = 8
#Less sensitive 
bowtie2_options = ['-D', '10', '-R', '1', '-N', '0', '-L', '20','-i', 'S,1,0.50', '-p', str(processors), '--reorder', '--no-hd', '--met', '20']
#very sensitive
#bowtie2_options = { '-D', '20', '-R', '3', '-N', '0', '-L', '20','-i', 'S,1,0.50','-p', num2str(params.processors)};
