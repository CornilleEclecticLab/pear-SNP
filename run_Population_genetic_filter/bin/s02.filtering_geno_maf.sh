#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s02.filtering_geno_maf
#SBATCH -o s02.filtering_geno_maf.out
#SBATCH -e s02.filtering_geno_maf.err
module purge
module load bioinfo/bcftools-1.14

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 \
--threads 8 \
-o $OUTPUT/pear.Combine_Chr.geno20_maf005.vcf.gz $OUTPUT/pear.Combine_Chr.vcf.gz

