#!/bin/bash
# Adaptation by Yuqi NIE
# 9/11/2022

#SBATCH -J s06.plink_vcf2bed
#SBATCH -o s06.plink_vcf2bed.NOGIT.%j.out
#SBATCH -e s06.plink_vcf2bed.NOGIT.%j.err
module purge
module load plink

source s00.config.sh

plink --vcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz \
--make-bed \
--const-fid \
--out $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.thin8k \
--set-missing-var-ids \@\:# \
--keep-allele-order \
--allow-extra-chr
