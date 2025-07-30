#!/usr/bin/env bash

#SBATCH -o s06.x1.log/s06.deleterious_mutation_sweep.02_merge.%J.out
#SBATCH -e s06.x1.log/s06.deleterious_mutation_sweep.02_merge.%J.err
#SBATCH --mem=4G

## This script would be submitted to slurm by s06.deleterious_mutation_sweep.01_pre.sh

set -euo pipefail

OUTDIR=../output/s06.dm_sweep

merge_burden() {
  local suffix=$1   # e.g. "CDS" or "synonymous" or "SNP"
  local prefix=$2   # e.g. "region.${suffix}.pop_mutation_burden.Rdata.txt"

  echo "Merging results for ${suffix}..."

  local combined="${OUTDIR}/${prefix}"
  # take header from any control file
  head -n1 "${OUTDIR}/control.${suffix}.pop_mutation_burden.txt" > "$combined"

  for f in "${OUTDIR}"/*.${suffix}.pop_mutation_burden.txt; do
    if [[ "$f" != "$combined" ]]; then
      tail -n +2 "$f" >> "$combined"
    fi
  done
}

merge_burden CDS         region.CDS.pop_mutation_burden.Rdata.txt
merge_burden synonymous  region.synonymous.pop_mutation_burden.Rdata.txt
merge_burden SNP         region.SNP.pop_mutation_burden.Rdata.txt

echo "All merges done at $(date)"
