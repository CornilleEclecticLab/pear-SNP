#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s05.filter_thin
#SBATCH -o s05.filter_thin.out
#SBATCH -e s05.filter_thin.err

module purge
module load bioinfo/samtools-1.14
module load bioinfo/vcftools-0.1.15

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

vcftools --gzvcf $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.vcf.gz --thin 8000 \
--recode --out $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k

mv $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.recode.vcf $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf

bgzip $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf &&\
tabix $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf.gz \
