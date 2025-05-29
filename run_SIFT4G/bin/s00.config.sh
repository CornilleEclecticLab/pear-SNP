# Set path, don't need to change
SCRIPT_DIR="$(realpath "$(dirname "${BASH_SOURCE[0]}")")"
PROJECT_ROOT="$(dirname ${SCRIPT_DIR})"
CONTAINER="${SCRIPT_DIR}/SIFT4G_Create_Genomic_DB.NOGIT./sift4g_db.sif"
INPUT_DIR="${PROJECT_ROOT}/input"
PROTEIN_DB="${PROJECT_ROOT}/database/uniref90.fasta"


############ Adjustable parameters ############

### for the first genome
# # Set input file names, adjust as needed
# GENOME_FASTA_FULL="${INPUT_DIR}/PyrusCommunis_BartlettDHv2.0.fasta"  # expected .fasta extension
# GFF_FILE="${INPUT_DIR}/PyrusCommunis_BartlettDHv2.0.gff"             # expected .gff extension, it will be transfrom into compressed gtf file 
# CHR_LIST="${INPUT_DIR}/PyrusCommunis_BartlettDHv2.0.chr_list.txt"    # One chromosome per line, no header, no blank lines

# # Set path, adjust as needed
# OUTPUT_DIR="${PROJECT_ROOT}/output/s01_database_PCOM"       # Adjust as needed
# CONFIG_FILE="${SCRIPT_DIR}/sift_config_PCOM.txt"            # Adjust as needed

# ORG="PyrusCommunis_BartlettDH"
# ORG_VERSION="Pcom_BartlettDHv20"


### for the second genome
# Set input file names, adjust as needed
GENOME_FASTA_FULL="${INPUT_DIR}/GWHBAOS00000000.genome.fasta"       # expected .fasta extension
GFF_FILE="${INPUT_DIR}/GWHBAOS00000000.no_blank.gff"                 # expected .gff extension, it will be transfrom into compressed gtf file 
CHR_LIST="${INPUT_DIR}/GWHBAOS00000000.chr_list.txt"

# Set path, adjust as needed
OUTPUT_DIR="${PROJECT_ROOT}/output/s01_database_PPY_long"       # Adjust as needed
CONFIG_FILE="${SCRIPT_DIR}/sift_config_PPY_long.txt"            # Adjust as needed

ORG="Pyrus_pyrifolia_Cuiguan"
ORG_VERSION="Ppy_Cuiguan"


############ END of Adjustable parameters ############


# Don't need to change
GENOME_FASTA=${GENOME_FASTA_FULL%.fasta}.chr.fa.gz
GTF_FILE="${GFF_FILE%.gff}.chr.gtf.gz"
