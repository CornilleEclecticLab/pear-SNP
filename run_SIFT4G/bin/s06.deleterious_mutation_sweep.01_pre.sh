#!/bin/bash

#SBATCH -J s06.deleterious_mutation_sweep.01_pre.sh
#SBATCH -o s06.log/s06.deleterious_mutation_sweep.01_pre.sh.%A_%a.%J.out
#SBATCH -e s06.log/s06.deleterious_mutation_sweep.01_pre.sh.%A_%a.%J.err
#SBATCH -c 1
#SBATCH --mem=48G
##SBATCH --array=1-18
#SBATCH --array=7-9


set -euo pipefail

module load conda
source activate intervaltree


# Common arguments
OUTDIR=../output/s06.dm_sweep
mkdir -p "$OUTDIR"


BEDS=()
for bed in ../input/selection_bed/selection.merge_across_methods.*.bed; do
  BEDS+=("$bed")
done


# Determine REGION and DEN based on task ID
case "${SLURM_ARRAY_TASK_ID}" in
  1|2|3) DEN=CDS ;;
  4|5|6) DEN=SNP ;;
  7|8|9) DEN=synonymous;;
  10|11|12) DEN=tolerated_SNP;;
  13|14|15) DEN=synonymous_allele;;
  16|17|18) DEN=tolerated_allele;;
  *)
    echo "Invalid task ID" >&2
    exit 1
    ;;
esac

case "${SLURM_ARRAY_TASK_ID}" in
  1|4|7|10|13|16) REGION=selection ;;
  2|5|8|11|14|17) REGION=control   ;;
  3|6|9|12|15|18) REGION=all_masked;;
esac

case "${SLURM_ARRAY_TASK_ID}" in
  1|2|3|4|5|6|7|8|9) SUMMARY="True";;
  *) SUMMARY="False";;
esac

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
    $SUMMARY \
    $BED_ARGS

echo "[$(date)] Finished ${REGION}/${DEN}"


# Submit the next job to merge results
if [[ "$SLURM_ARRAY_TASK_ID" -eq 7 ]]; then
  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} s06.x2.deleterious_mutation_sweep.02_merge.sh
fi

echo "[$(date)] Running ${REGION}/${DEN} on task ${SLURM_ARRAY_TASK_ID}"
