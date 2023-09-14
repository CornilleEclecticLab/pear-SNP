#!/bin/bash
#SBATCH -J s01.get_pca_plink
#SBATCH -o s01.get_pca_plink.%J.out
#SBATCH -e s01.get_pca_plink.%J.err
module purge
module load plink/1.90b6.18

WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_PCA_plink/"
set="pear_Jul2023_noclone.Combine_Chr.geno20_maf005.anno.syno.thin8k"
OUTDIR="$WORKDIR/output/branch8.batch15.pear_Jul2023_noclone"
INPUT="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/branch8.batch15.pear_Jul2023_noclone/pear_Jul2023_noclone.Combine_Chr.geno20_maf005.anno.syno.thin8k"

mkdir -p $OUTDIR


plink -bfile $INPUT \
    --allow-extra-chr \
    --pca var-wts \
    --make-rel \
    --out ${set}.pca


# perl get_R_data.pl # Disused, as Plink couldn't calculate the Proportion of variance explained(pve) directly, but this script set eigenvalues as PVE. To calculate PVE, use following new script.
perl bin.calculate_pve.pl $set


#mv ./plink.* $OUTDIR
#mv ./*.eigenv* $OUTDIR
#mv ./pca_input.data $OUTDIR
#mv ./${set}.pac* $OUTDIR
