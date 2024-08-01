#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02

# Update on 2024-05-22 by Yuqi:
# 1. change to use plink to thin the vcf file rather vcftools, because vcftools would change the format of the vcf file, e.g missing values from "./." to "."

#SBATCH -J s05.filter_thin
#SBATCH -o s05.filter_thin.NOGIT.%j.out
#SBATCH -e s05.filter_thin.NOGIT.%j.err

module purge
module load samtools/1.14
module load plink/1.90b6.18
source s00.config.sh
# source ~/.bashrc
# conda activate vcftools-0.1.15



# vcftools --gzvcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf.gz --thin 8000 \
#     --recode --out $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k

# mv $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.recode.vcf \
#     $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf

# bgzip $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf &&\
# tabix $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz 


plink --vcf "$OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf.gz" \
--make-bed \
--recode vcf-iid \
--const-fid \
--bp-space 8000 \
--out $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k \
--set-missing-var-ids \@\:# \
--keep-allele-order \
--allow-extra-chr

bgzip -k $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf
tabix -p vcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz
