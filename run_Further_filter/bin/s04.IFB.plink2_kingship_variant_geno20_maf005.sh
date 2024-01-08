#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship_variant_geno20.%J.out
#SBATCH -e s04.kinship_variant_geno20.%J.err
#SBATCH -c 8

module purge
module load plink2/2.00a2.3
module load plink/1.90b6.18
module load bcftools/1.14


BATCH='pear_Oct2023'


BIN=$(pwd)
WORKDIR=$(dirname $BIN)
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/s04.${BATCH}_variant_geno20_maf005"
PREFIX="$BATCH.Combine_chr_variant_geno20_maf005"
OUTPUTs04=$OUTPUT/$PREFIX.s04.clone_filtering

mkdir -p $OUTPUT
mkdir -p $OUTPUTs04

bcftools concat -f $INPUT/s04.${BATCH}.chr.variant.vcf.list \
    | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' \
    -O z4 \
    --threads 8 \
    -o $OUTPUT/$PREFIX.vcf.gz


plink --vcf $OUTPUT/$PREFIX.vcf.gz \
--make-bed \
--const-fid \
--out $OUTPUT/$PREFIX \
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr 



# Do not change --bfile to --vcf 
    # If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
    # Invalid chromosome code 'Pp01' on line 240 of --vcf file.
    # (Use --allow-extra-chr to force it to be accepted.)
plink2 --bfile $OUTPUT/$PREFIX \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out $OUTPUTs04/without_clone0.354 


plink2 --bfile $OUTPUT/$PREFIX \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.177 \
--make-bed \
--out $OUTPUTs04/without_1st-degree_relation0.177 


plink2 --bfile $OUTPUT/$PREFIX \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.088 \
--make-bed \
--out $OUTPUTs04/without_2nd-degree_relation0.088 



