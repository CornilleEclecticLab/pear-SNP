#!/bin/bash
#SBATCH -J s01.get_pca_plink
#SBATCH -o s01.get_pca_plink.out
#SBATCH -e s01.get_pca_plink.err
module purge
module load bioinfo/plink-v1.90b5.3


WORKDIR="/home/zruilin/work/Peach/run_PCA_plink"
INPUT="$WORKDIR/input/Combine_Chr.geno20_maf005.anno.syno.thin8k"
OUTDIR="$WORKDIR/output/s01.get_pca_plink"
mkdir $OUTDIR


plink -bfile $INPUT --allow-extra-chr --pca

perl get_R_data.pl

mv ./*.eigenv* $OUTDIR/
mv ./pca_input.data $OUTDIR/

