#!/bin/bash
# by Yuqi NIE at 2023-09-15 16:10:31

# Acknowledgement:
# 	Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861


#SBATCH -J s03.extract_pure_sample_filter_geno_maf_to_nex
#SBATCH -o s03.extract_pure_sample_filter_geno_maf_to_nex.%j.out
#SBATCH -e s03.extract_pure_sample_filter_geno_maf_to_nex.%j.err
module purge
module load bcftools/1.14


# Define input and output paths and filenames
# WORK_DIR="/shared/home/ynie/work/pear/run_SVDQuartets" # Replace with the actual path to your input directory
WORK_DIR="../"
PREFIX="pear_Jul2023_noclone"   # Replace with your desired prefix
MIDDLE="Combine_Chr.geno20_maf005.anno.syno.thin8k"
INPUT_VCF="$WORK_DIR/input/$PREFIX.$MIDDLE.vcf.gz"
OUTPUT_DIR="$WORK_DIR/output/s03.extract_pure_sample_filter_geno_maf_to_nex"
OUTPUT_VCF="$OUTPUT_DIR/$PREFIX.$MIDDLE.Nomix.vcf.gz"
S02_OUTPUT_DIR="$WORK_DIR/output/s02.choose_pure_samples_manually"

# Define the sample list file
SAMPLE_LIST="$S02_OUTPUT_DIR/s02.non_admixed_list.txt"  # or SAMPLE_LIST="^PATH", where ^ means exclude
RENAME_LIST="$S02_OUTPUT_DIR/s02.non_admixed_rename_list.txt"  # or SAMPLE_LIST="^PATH", where ^ means exclude

# Define filtering criteria
FILTER_EXPRESSION='F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN'

# Define the number of threads for parallel processing
THREADS=8

echo 'Your input file is: ' $INPUT_VCF

# Make sure the folder exist
mkdir -p $OUTPUT_DIR 

# Use bcftools to filter and process the VCF file
bcftools view "$INPUT_VCF" \
    --samples-file "$SAMPLE_LIST" \
| bcftools reheader \
    --samples "$RENAME_LIST" \
| bcftools filter \
    --exclude "$FILTER_EXPRESSION" \
    --output-type z4 \
    --threads "$THREADS" \
    --output "$OUTPUT_VCF"


VCF2PHYLIP="vcf2phylip.py"

if [ ! -e "$VCF2PHYLIP" ]; then
	wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
fi


# Convert VCF to NEX
/usr/bin/python3 vcf2phylip.py \
	--input $OUTPUT_VCF \
    --output-folder $OUTPUT_DIR \
    --phylip-disable \
	--nexus


# Combine nexus files
cat $OUTPUT_DIR/$PREFIX.$MIDDLE.Nomix.*.nexus "$S02_OUTPUT_DIR/s02.taxpartitions.nex" \
    > $OUTPUT_DIR/$PREFIX.$MIDDLE.Nomix.parts.nex

echo 'Your output file is: ' 
echo  "$OUTPUT_DIR/$PREFIX.$MIDDLE.Nomix.parts.nex"
