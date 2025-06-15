#!/bin/bash


#SBATCH -J s06.deleterious_mutation_sweep.01_pre.sh
#SBATCH -o s06.deleterious_mutation_sweep.01_pre.sh.%A_%a.out
#SBATCH -e s06.deleterious_mutation_sweep.01_pre.sh.%A_%a.err
#SBATCH -c 1
#SBATCH --mem=48G
#SBATCH --array=1-9


set -euo pipefail

module load conda
source activate intervaltree


# Common arguments
OUTDIR=../output/s06.dm_sweep
mkdir -p "$OUTDIR"

BEDS=(
  ../input/selection_bed/selection.outliers.betu.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.cauc.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.comm_Dessert.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.comm_Perry.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.pyra.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.pyri_JP.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.Sand_CN-SE.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.Sand_CN-SW.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.ussu.filtered.sorted.merged.bed
  ../input/selection_bed/selection.outliers.White.filtered.sorted.merged.bed
)

# Determine REGION and DEN based on task ID
case "${SLURM_ARRAY_TASK_ID}" in
  1|2|3) DEN=CDS ;;
  4|5|6) DEN=SNP ;;
  7|8|9) DEN=synonymous ;;
  *) 
    echo "Invalid task ID" >&2
    exit 1
    ;;
esac

case "${SLURM_ARRAY_TASK_ID}" in
  1|4|7) REGION=selection ;;
  2|5|8) REGION=control   ;;
  3|6|9) REGION=all_masked;;
esac

# Submit the next job to merge results
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} s06.deleterious_mutation_sweep.02_merge.sh
fi

echo "[$(date)] Running ${REGION}/${DEN} on task ${SLURM_ARRAY_TASK_ID}"

# Build BED_ARGS if needed
if [[ "$REGION" == "all_masked" ]]; then
  BED_ARGS=""
else
  BED_ARGS="${BEDS[*]}"
fi

# Run the Python script
python3 s06.x1.dm_extract_sweep.py \
    "$OUTDIR" \
    "$REGION" \
    "$DEN" \
    $BED_ARGS

echo "[$(date)] Finished ${REGION}/${DEN}"
