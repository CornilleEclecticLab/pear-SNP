# This config file is for the pixy on the IFB cluster.
LOAD_PIXY="module load pixy/1.2.7.beta1"
LOAD_BCFTOOLS="module load bcftools/1.14"

# The following sets only effect in the s01.run_pixy.sh script for testing,
# not in s02.run_pixy_by_chromosome.py
INPUT_VCF="test.vcf.gz"
BIN_DIR=$(pwd)
INPUT_DIR="${BIN_DIR}/../input"
OUTPUT_DIR="${BIN_DIR}/../output"

PREFIX="Ag1000_sampleIDs_popfile"
INPUT_POP="${INPUT_DIR}/${PREFIX}.txt"
