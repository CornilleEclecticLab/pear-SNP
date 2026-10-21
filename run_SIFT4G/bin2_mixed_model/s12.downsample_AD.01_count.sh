#!/bin/bash

#SBATCH -J s12.downsample_AD.01_count
#SBATCH -o s12.log/s12.downsample_AD.01_count.%A_%a.out
#SBATCH -e s12.log/s12.downsample_AD.01_count.%A_%a.err
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --array=1-34

# In silico downsampling of allelic read counts to a common depth.

# Usage: sbatch s12.downsample_AD.01_count.sh
# Then:  python3 s12.x2.merge_counts.py

set -euo pipefail

source s00.load_bcftools.sh
source s00.load_bedtools.sh

TARGET=10
REPS=3

OUTDIR=../output/s12.downsample_AD
TMPDIR_S10=../temp/s12.downsample_AD
mkdir -p "$OUTDIR/counts" "$TMPDIR_S10"

# Task list: PPY, PCOM
TASK_LINE=$(for ref in PPY PCOM; do
              awk -v r=$ref '{print r"\t"$1}' ../input/s03.variant_vcf_list.$ref.txt
            done | sed -n "${SLURM_ARRAY_TASK_ID}p")
REF=$(cut -f1 <<< "$TASK_LINE")
VCF=$(cut -f2 <<< "$TASK_LINE")
CHR=$(basename "$VCF" | cut -d. -f2)
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: $REF $CHR $VCF"

# Mapped depth
DEPTH=$TMPDIR_S10/$REF.depth.txt
if [[ "$REF" == "PPY" ]]; then
  awk -F'\t' 'NR>1{print $1"\t"$4}' depth.txt > $DEPTH.$CHR
else
  awk -F'\t' 'NR>1{print $2"\t"$5}' \
    ../../run_GATK_variant_calling_refer_on_communis/output/s02.stat_format.pear.ref_comm.Jul2024.txt > $DEPTH.$CHR
fi

BED=$TMPDIR_S10/$REF.$CHR.cds_masked.merged.bed
awk -v c=$CHR '$1==c' ../input/$REF.cds_masked.bed | sort -k2,2n | bedtools merge > $BED

SAMPLES=$TMPDIR_S10/$REF.$CHR.samples.txt
cut -f1 ../input/s01.sample_tab_population_$REF.txt > $SAMPLES.in
bcftools view -h -s $(paste -sd, $SAMPLES.in) "$VCF" | bcftools query -l > $SAMPLES

bcftools view -s $(paste -sd, $SAMPLES) -R $BED "$VCF" -Ou \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%AD:%DP]\n' \
| /shared/projects/pear_snp2/conda/env/pixy2/bin/python3 s12.x1.downsample_AD_count.py \
    --chrom $CHR \
    --samples $SAMPLES \
    --depth $DEPTH.$CHR \
    --sift ../output/s02_annot/$REF/merged.SIFTannotations.xls.cleaned \
    --target $TARGET \
    --reps $REPS \
    --out $OUTDIR/counts/$REF.$CHR.counts.tsv

rm -f $DEPTH.$CHR $SAMPLES.in
echo "[$(date)] Finished $REF $CHR"
