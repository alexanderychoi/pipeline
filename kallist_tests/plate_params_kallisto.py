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

####Path to fastq files in the given directory:
#dir_path_fastqs = '../fastqs_sisi/fastqs'
#dir_path_alignment = '../fastqs_sisi/alignment/'
dir_path_fastqs = '../fastqs_kallisto/fastqs/'
dir_path_alignment = '../fastqs_kallisto/alignment/'

####Path to the fasta file [Kallisto]
reference_genome = '../reference_genomes/mm9/Transcriptome/transcripts.fa'

####Path to the .idx file (will be created if does not exist)
reference_idx = '../reference_genomes/mm9/Transcriptome/transcripts.idx'

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
umi_length = 5

####TAC length:
tac_length = 3

####Number of G's to delete
g_length = 4

####Tso sequence
tso = 'AAGCAGTGGTATCAACGCAGAGTAC'

####Number of Ts or As in a row beyond which the polyT or polyA sequence is trimmed
num_ta = 10
error_string_t = 'T'*num_ta
error_string_a = 'A'*num_ta

####Minimum length for RNA seq (if < minimum_length,  the read is dismisse)
minimum_length = 30


#######################
##   Kallisto vars   ##
#######################

####

