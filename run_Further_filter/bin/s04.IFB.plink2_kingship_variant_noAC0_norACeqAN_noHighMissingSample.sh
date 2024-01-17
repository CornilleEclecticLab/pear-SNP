#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship.NOGIT.%J.out
#SBATCH -e s04.kinship.NOGIT.%J.err
#SBATCH -c 8

module purge
module load plink2/2.00a2.3
module load bcftools/1.14


# The $BATCH is revelent to the input file of vcf list named
# "../input/s04.${BATCH}.chr.variant.vcf.list"
BATCH='pear_Dec2023'


BIN=$(pwd)
WORKDIR=$(dirname $BIN)
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/s04.${BATCH}.variant_noAC0_norACeqAN_noHighMissingSample"
PREFIX="$BATCH.Combine_chr"
OUTPUTs04=${OUTPUT}/${PREFIX}.s04.clone_filtering

mkdir -p ${OUTPUT}
mkdir -p ${OUTPUTs04}

bcftools concat -f ${INPUT}/s04.${BATCH}.chr.variant.vcf.list \
    --threads 8 \
    | bcftools view \
        -S ^${WORKDIR}/output/s04.pear_Dec2023_variant_noAC0_norACeqAN/pear_Dec2023.Combine_chr.vcf.gz.stats.txt.higher40_missing_rate_sample.list.txt \
        --threads 8 \
    | bcftools filter -e 'AC==0 || AC==AN' \
        -O z4 \
        --threads 8 \
        -o ${OUTPUT}/${PREFIX}.vcf.gz

bcftools stats -s - ${OUTPUT}/${PREFIX}.vcf.gz \
    > ${OUTPUT}/${PREFIX}.vcf.gz.stats.txt

cat ${OUTPUT}/${PREFIX}.vcf.gz.stats.txt \
    | grep -E '^PSC' \
    | cut -f 3,14 \
    > ${OUTPUT}/${PREFIX}.vcf.gz.ind_miss.txt

cat ${OUTPUT}/${PREFIX}.vcf.gz.stats.txt \
    | grep -v '#' \
    | grep 'number of SNPs:' \
    | cut -f 3,4 \
    > ${OUTPUT}/${PREFIX}.vcf.gz.SNPs_number.txt

# Now we don't use plink to transform vcf to bed format.
# plink --vcf $OUTPUT/$PREFIX.vcf.gz \
# 	--make-bed \
# 	--const-fid \
# 	--out $OUTPUT/$PREFIX \
# 	--set-missing-var-ids @:# \
# 	--keep-allele-order \
# 	--allow-extra-chr 


# Invalid Annotation:
    # Do not change --bfile to --vcf 
    # If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
    # Invalid chromosome code 'Pp01' on line 240 of --vcf file.
    # (Use --allow-extra-chr to force it to be accepted.)
plink2 --vcf $OUTPUT/$PREFIX.vcf.gz \
    --allow-extra-chr \
    --make-king-table \
    --king-cutoff 0.354 \
    --make-bed \
    --out $OUTPUTs04/without_clone0.354

plink2 --vcf $OUTPUT/$PREFIX.vcf.gz \
    --allow-extra-chr \
    --make-king-table \
    --king-cutoff 0.177 \
    --make-bed \
    --out $OUTPUTs04/without_1st-degree_relation0.177

plink2 --vcf $OUTPUT/$PREFIX.vcf.gz \
    --allow-extra-chr \
    --make-king-table \
    --king-cutoff 0.088 \
    --make-bed \
    --out $OUTPUTs04/without_2nd-degree_relation0.088
