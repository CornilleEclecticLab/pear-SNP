#!/bin/bash
#SBATCH -J s01.get_pca_plink
#SBATCH -o s01.get_pca_plink.%J.out
#SBATCH -e s01.get_pca_plink.%J.err
module purge
module load plink/1.90b6.18

source s00.config.sh

mkdir -p "$OUTDIR"


plink -bfile "${INPUT}" \
    --allow-extra-chr \
    --pca var-wts \
    --make-rel \
    --out "${OUTDIR}/${set}.pca"

cd "$OUTDIR" || exit 1
# perl get_R_data.pl # Disused, as Plink couldn't calculate the Proportion of variance explained(pve) directly, but this script set eigenvalues as PVE. To calculate PVE, use following new script.
perl "${WORKDIR}/bin/bin.calculate_pve.pl" "$set"


#mv ./plink.* $OUTDIR
#mv ./*.eigenv* $OUTDIR
#mv ./pca_input.data $OUTDIR
#mv ./${set}.pac* $OUTDIR
