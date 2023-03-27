#!/bin/bash
# Adaptation by Yuqi NIE
# 9/11/2022

#SBATCH -J s06.plink_vcf2bed
#SBATCH -o s06.plink_vcf2bed.out
#SBATCH -e s06.plink_vcf2bed.err
module purge
module load plink

WORKDIR="/shared/ifbstor1/projects/pear_snp2/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.whole_pear"
PREFIX="whole_pear"

plink --vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf.gz \
--make-bed \
--const-fid \
--out $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.thin8k \
--set-missing-var-ids \@\:# \
--keep-allele-order \
--allow-extra-chr
