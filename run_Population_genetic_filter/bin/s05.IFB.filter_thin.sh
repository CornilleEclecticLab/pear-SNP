#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

#SBATCH -J s05.filter_thin
#SBATCH -o s05.filter_thin.NOGIT.%j.out
#SBATCH -e s05.filter_thin.NOGIT.%j.err

module purge
module load samtools/1.14
source s00.config.sh
source ~/.bashrc
conda activate vcftools-0.1.15

vcftools --gzvcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf.gz --thin 8000 \
    --recode --out $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k

mv $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.recode.vcf \
    $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf

bgzip $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf &&\
tabix $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz 
