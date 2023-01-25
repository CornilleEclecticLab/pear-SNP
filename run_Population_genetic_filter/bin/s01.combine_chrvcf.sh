#!/bin/bash
# By Xilong CHEN
# Create Date: 8/21/2020 12:19
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2022-01-23

#SBATCH -J s01.combine_chrvcf
#SBATCH -o s01.combine_chrvcf.out
#SBATCH -e s01.combine_chrvcf.err

module purge
module load bioinfo/bcftools-1.14

WORKDIR="/work/zruilin/Peach/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

bcftools concat -f $INPUT/chr_variant_file_list.txt \
-O z \
-o $OUTPUT/Combine_Chr.vcf.gz
