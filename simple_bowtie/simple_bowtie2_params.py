####Path to .fastq files
dir_path_fastqs = '../data/070215_drop/fastqs/'

####Path to .sam files
dir_path_alignment = '../data/070215_drop/alignment/'

####Path to the reference genome files
#reference_genome = '../reference_genomes/mm9/Transcriptome/transcriptome'
reference_genome = '../reference_genomes/Human_Mouse/transcriptome'


#######################
## Preprocessing var ##
#######################

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
bowtie2_options = ['-D', '10', '-R', '1', '-N', '0', '-L', '20','-i', 'S,1,0.50', '--local', '-p', str(processors), '--no-hd', '--met', '20']
#very sensitive
#bowtie2_options = { '-D', '20', '-R', '3', '-N', '0', '-L', '20','-i', 'S,1,0.50','-p', num2str(params.processors)};