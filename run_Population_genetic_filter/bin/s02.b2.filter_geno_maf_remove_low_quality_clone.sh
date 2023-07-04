#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s02.b2.filtering_geno_maf_remove_low_clone
#SBATCH -o s02.b2.filtering_geno_maf_remove_low_clone.%j.out
#SBATCH -e s02.b2.filtering_geno_maf_remove_low_clone.%j.err

module purge
module load bioinfo/bcftools-1.14

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.whole_pear"
PREFIX="whole_pear"

bcftools view \
	-S ^s02.b2.remove_low_quality_clone.list.txt \
	$OUTPUT/$PREFIX.Combine_Chr.vcf.gz \
|bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 \
	-o $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz # --threads 8 

