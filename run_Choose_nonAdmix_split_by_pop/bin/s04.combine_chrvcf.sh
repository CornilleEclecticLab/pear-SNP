#!/bin/bash
# By Xilong CHEN
# Create Date: 8/21/2020 12:19
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2023-11-29 16:09:06

#SBATCH -J s04.combine_chrvcf
#SBATCH -o s04.combine_chrvcf.NOGIT.%J.out
#SBATCH -e s04.combine_chrvcf.NOGIT.%J.err
#SBATCH -c 4

module purge
module load bcftools/1.14

BATCH="pear_July2023"

BINDIR=$(pwd)
WORKDIR=$(dirname $BINDIR)
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/s04.combine_chrvcf"
OUTPUT_PREFIX="filtered_pixy_variant.non_admix"

mkdir -p $OUTPUT

bcftools concat -f ${INPUT}/s04.chr_vcf_list.txt \
    -O z4 \
    --threads 4 \
    -o ${OUTPUT}/${BATCH}.${OUTPUT_PREFIX}.vcf.gz

tabix -p vcf ${OUTPUT}/${BATCH}.${OUTPUT_PREFIX}.vcf.gz

bcftools stats ${OUTPUT}/${BATCH}.${OUTPUT_PREFIX}.vcf.gz \
> ${OUTPUT}/${BATCH}.${OUTPUT_PREFIX}.vcf.gz.stats
