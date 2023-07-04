#!/bin/bash
# By Xilong CHEN
# Create Date: 8/21/2020 12:19
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2022-01-23

#SBATCH -J s01.combine_chrvcf
#SBATCH -o s01.combine_chrvcf.out
#SBATCH -e s01.combine_chrvcf.err

module purge
module load bcftools/1.14

WORKDIR="/shared/ifbstor1/projects/pear_snp2/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

bcftools concat -f $INPUT/whole_pear.cvr.variant.vcf.list \
-O z4 \
--threads 8 \
-o $OUTPUT/branch7.whole_pear/whole_pear.Combine_Chr.vcf.gz
