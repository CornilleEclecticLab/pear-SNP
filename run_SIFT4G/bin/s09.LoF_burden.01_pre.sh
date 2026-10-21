#!/bin/bash

#SBATCH -J s09.LoF_burden.01_pre.sh
#SBATCH -o s09.log/s09.LoF_burden.01_pre.sh.%A_%a.out
#SBATCH -e s09.log/s09.LoF_burden.01_pre.sh.%A_%a.err
#SBATCH -c 1
#SBATCH --mem=48G
#SBATCH --array=1-9

#Description:
    # LoF burden (denominator = synonymous SNP) per region and LoF set.

    # Requires s08.calculate_LoF.sh and s09.a1.extract_synonymous.sh done for PCOM and PPY.

set -euo pipefail

for REF in PCOM PPY; do
  for f in merged.LoF_merged.tsv merged.synonymous_consensus.tsv; do
    if [[ ! -s ../output/s08.calculate_LoF/${REF}/${f} ]]; then
      echo "Missing ../output/s08.calculate_LoF/${REF}/${f}: run s08.calculate_LoF.sh / s09.a1.extract_synonymous.sh ${REF} first" >&2
      exit 1
    fi
  done
done

module load conda
source activate intervaltree


BEDS=()
for bed in ../input/selection_bed/selection.merge_across_methods.*.bed; do
  BEDS+=("$bed")
done

case "${SLURM_ARRAY_TASK_ID}" in
  1|2|3) LOF_SET=consensus ;;
  4|5|6) LOF_SET=union ;;
  7|8|9) LOF_SET=hc ;;
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

if [[ "$REGION" == "all_masked" ]]; then
  BED_ARGS=""
else
  BED_ARGS="${BEDS[*]}"
fi

OUTDIR=../output/s09.LoF_burden/${LOF_SET}
mkdir -p "$OUTDIR"

echo "[$(date)] Running ${LOF_SET}/${REGION} on task ${SLURM_ARRAY_TASK_ID}"

python3 s09.x1.LoF_burden.py \
    "$OUTDIR" \
    "$REGION" \
    "$LOF_SET" \
    $BED_ARGS

echo "[$(date)] Finished ${LOF_SET}/${REGION}"


# Submit the next job
if [[ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MIN" ]]; then
  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} s09.LoF_burden.02_merge.sh
fi
