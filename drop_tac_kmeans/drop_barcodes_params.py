######################################
## drop_barcodes_params.py allows you to change
## the parameters for creating the matrix of plasmid
## barcode reads


#### directory containing pickle files that include groups variable
#dir_path_pickle = '/home/iamcam/Documents/RNAseqdata/DS35/Brn-barcodes/pickle/'
#dir_path_pickle = '/home/iamcam/Documents/RNAseqdata/DS31/Pooled_barcodes/pickle/'
dir_path_pickle = '/home/iamcam/Documents/RNAseqdata/DS31/Sep_barcodes/pickle/'

#### directory containing barcode file directory
#dir_path_barcodes = '/home/iamcam/Documents/RNAseqdata/DS35/Brn-barcodes/'
#dir_path_barcodes = '/home/iamcam/Documents/RNAseqdata/DS31/Pooled_barcodes/'
dir_path_barcodes = '/home/iamcam/Documents/RNAseqdata/DS31/Sep_barcodes/'

#### Name of data matrix to get ordering of barcodes from 
#datafile = '/home/iamcam/Documents/RNAseqdata/DS35/brn2_matrix.txt'
#datafile = '/home/iamcam/Documents/RNAseqdata/DS31/PoolD2_matrix.txt'
datafile = '/home/iamcam/Documents/RNAseqdata/DS31/SeparateD2_matrix.txt'


#### String of 8 bases after plasmid barcodes
p_str = 'GCCTAGGG'

####Path to the reference genome files
reference_genome = '/home/iamcam/Documents/reference_genomes/barcodes/barcodes'

#### Path to the species fasta file
fasta_path = '/home/iamcam/Documents/reference_genomes/barcodes/barcodes.fa'
