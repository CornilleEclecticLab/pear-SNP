#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s05.filter_thin
#SBATCH -o s05.filter_thin.%j.out
#SBATCH -e s05.filter_thin.%j.err

module purge
module load samtools/1.14
. ~/.bashrc

WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch8.batch15.pear_Jul2023_noclone"
PREFIX="pear_Jul2023_noclone"

conda activate vcftools-0.1.15
vcftools --gzvcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf.gz --thin 8000 \
    --recode --out $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k

mv $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.recode.vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf

conda activate tabix-0.2.5
bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf &&\
tabix $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf.gz \
