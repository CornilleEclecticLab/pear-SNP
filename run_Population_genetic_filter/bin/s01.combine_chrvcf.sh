#!/bin/bash
# By Xilong CHEN
# Create Date: 8/21/2020 12:19
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2022-01-23

#SBATCH -J s01.combine_chrvcf
#SBATCH -o s01.combine_chrvcf.NOGIT.%J.out
#SBATCH -e s01.combine_chrvcf.NOGIT.%J.err
#SBATCH -c 4

module purge
module load bcftools/1.14

WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch15.pear_Oct2023"

mkdir -p $OUTPUT

bcftools concat -f $INPUT/pear_Oct2023.chr.variant.vcf.list \
-O z4 \
--threads 8 \
-o $OUTPUT/pear_Oct2023.Combine_Chr.vcf.gz
