#!/usr/bin/env bash

#SBATCH -o s09.log/s09.LoF_burden.02_merge.%J.out
#SBATCH -e s09.log/s09.LoF_burden.02_merge.%J.err
#SBATCH --mem=4G

## This script would be submitted to slurm by s09.LoF_burden.01_pre.sh

set -euo pipefail

merge_regions() {
  local outdir=$1   # e.g. ../output/s09.LoF_burden/consensus
  local suffix=$2   # e.g. "pop_mutation_burden" or "pop_mutation_burden.summary"

  local combined="${outdir}/region.synonymous.${suffix}.Rdata.txt"
  head -n1 "${outdir}/control.synonymous.${suffix}.txt" > "$combined"

  for region in all_masked control selection; do
    tail -n +2 "${outdir}/${region}.synonymous.${suffix}.txt" >> "$combined"
  done
}

for LOF_SET in consensus union hc; do
  [[ -d ../output/s09.LoF_burden/${LOF_SET} ]] || continue
  echo "Merging results for ${LOF_SET}..."
  merge_regions ../output/s09.LoF_burden/${LOF_SET} pop_mutation_burden
  merge_regions ../output/s09.LoF_burden/${LOF_SET} pop_mutation_burden.summary
done

echo "All merges done at $(date)"
