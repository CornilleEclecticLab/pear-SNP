#!/bin/bash
# Adaptation based on the script by Xilong CHEN
# Edited by Yuqi NIE at 2023-09-27 13:39:03

#SBATCH -J s02.filtering
#SBATCH -o s02.IFB.filter_SNP_for_PoPLDdecay.%J.out
#SBATCH -e s02.IFB.filter_SNP_for_PoPLDdecay.%J.err
module purge
module load bcftools/1.16

WORKDIR="$(pwd)/.."
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch8.pear_Jul2023"
PREFIX="pear_Jul2023"
BOUTPUT="$WORKDIR/output/branch8.batch15.pear_Jul2023_LDdecay"
BPREFIX="pear_Jul2023_LDdecay"

mkdir -p $BOUTPUT

# -S ^list
bcftools view \
    -S $INPUT/s02.keep_sample_list_for_PopLDdecay.txt \
    $INPUT/$PREFIX.Combine_Chr.vcf.gz \
| bcftools annotate \
    -x INFO,INFO,^FORMAT/GT \
| bcftools filter \
    -e 'F_MISSING > 0.25 || MAF <= 0.005' \
    -O z4 \
    --threads 8 \
    -o $BOUTPUT/$BPREFIX.Combine_Chr.geno25_maf0005_for_PopLDdecay.vcf.gz
