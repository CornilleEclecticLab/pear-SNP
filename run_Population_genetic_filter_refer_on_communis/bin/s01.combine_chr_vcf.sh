#!/bin/bash
# By Xilong CHEN
# Create Date: 8/21/2020 12:19
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2022-01-23
# Edit for running at 2024-08-01 15:41:36

#SBATCH -J s01.combine_chrvcf
#SBATCH -o s01.combine_chrvcf.NOGIT.%J.out
#SBATCH -e s01.combine_chrvcf.NOGIT.%J.err
#SBATCH -c 4

module purge
module load bcftools/1.14

WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter_refer_on_communis"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/pear_Jul2024_ref_comm.branch17.whole_pear"

mkdir -p $OUTPUT

bcftools concat -f $INPUT/pear.Jul2024.ref_on_communis.chr.variant.vcf.list.txt \
-O z4 \
--threads 8 \
-o $OUTPUT/pear_Jul2024_ref_comm.Combine_chr.vcf.gz

tabix -p vcf $OUTPUT/pear_Jul2024_ref_comm.Combine_chr.vcf.gz

