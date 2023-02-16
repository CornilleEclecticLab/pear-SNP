#!/bin/bash

#SBATCH -J s01.get_bed_plink
#SBATCH -o s01.get_bed_plink.out
#SBATCH -e s01.get_bed_plink.err
module purge
module load bioinfo/plink-v1.90b5.3

WORKDIR="/work/xchen/XCHEN_apple_genome/analysis_2021_Feb/Population_genetics_v2"

mkdir $WORKDIR/output/s03_fastSt_subset
INPUT="$WORKDIR/output/s02_keep_subset/Combine_Chr_rm_ind_anno_geno20_maf005_syno_LDpruned_thin5k_subset.vcf.gz"
OUTPUT="$WORKDIR/output/s03_fastSt_subset/Combine_Chr_rm_ind_anno_geno20_maf005_syno_LDpruned_thin5k_subset"

plink --vcf $INPUT \
--make-bed \
--const-fid \
--out $OUTPUT \
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr
