#!/bin/bash
#SBATCH -J s01.get_pca_plink
#SBATCH -o s01.get_pca_plink.out
#SBATCH -e s01.get_pca_plink.err
module purge
module load bioinfo/plink-v1.90b5.3


WORKDIR="/work/zruilin/pear/run_PCA_plink/"
INPUT="/work/zruilin/pear/run_Population_genetic_filter/output/batch1.reduce3EupropeWildPear/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k"
OUTDIR="$WORKDIR/output/batch1.reduce3EupropeWildPear"
set="pear.Combine_Chr.geno20_maf005.anno.syno.thin8k"

mkdir -p $OUTDIR


plink -bfile $INPUT \
    --allow-extra-chr \
    --pca var-wts \
    --make-rel \
    --out ${set}.pca

# perl get_R_data.pl
perl bin.calculate_pve.pl $set


#mv ./plink.* $OUTDIR
#mv ./*.eigenv* $OUTDIR
#mv ./pca_input.data $OUTDIR
#mv ./${set}.pac* $OUTDIR
