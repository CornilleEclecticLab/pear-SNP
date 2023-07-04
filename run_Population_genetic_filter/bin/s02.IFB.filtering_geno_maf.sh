#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s02.filtering_geno_maf
#SBATCH -o s02.filtering_geno_maf.%j.out
#SBATCH -e s02.filtering_geno_maf.%j.err
module purge
module load bcftools/1.14

WORKDIR="/shared/ifbstor1/projects/pear_snp2/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.whole_pear"
PREFIX="whole_pear"

bcftools view $OUTPUT/$PREFIX.Combine_Chr.vcf.gz \
	-S ^s02.b7.whole_pear.remove_low_quality_clone.list.txt \
|bcftools filter \
	-e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' \
	-O z4 \
	--threads 8 \
	-o $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz 

