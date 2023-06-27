#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s05.filter_thin
#SBATCH -o s05.filter_thin.%j.out
#SBATCH -e s05.filter_thin.%j.err

module purge
module load bioinfo/samtools-1.14
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.batch12.kinship_base_on_branch7_batch11_remove_clone_low_quality"
PREFIX="whole_pear"

vcftools --gzvcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf.gz --thin 8000 \
--recode --out $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k

mv $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.recode.vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf

bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf &&\
tabix $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf.gz \
