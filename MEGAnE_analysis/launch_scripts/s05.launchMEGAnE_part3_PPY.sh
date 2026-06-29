#!/bin/bash
#SBATCH --job-name=PPY_MEGAnE_part3
#SBATCH --output=../logs/PPY_MEGAnE_part3.out
#SBATCH --error=../logs/PPY_MEGAnE_part3.err
#SBATCH --time=6-00:00:00
#SBATCH -p compute
#SBATCH -c 32
#SBATCH --mem=64G

set -euo pipefail

module load singularity

sif=/scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif

# Base directory where per-sample outputs were written in Step 1
OUT_BASE="../output/PPY/MEGAnE_result"

# Species/cohort name variable
SPECIES='PPY'

################################################################################
#### Step 2. Joint calling

# 1) List all per-sample MEGAnE output directories (one level under OUT_BASE)
#    Write ABSOLUTE paths into dirlist.txt (some tools prefer absolute paths).
#    If you want to be stricter, you can filter only directories containing expected result files.
if [[ ! -d "$OUT_BASE" ]]; then
  echo "No output base directory found at $OUT_BASE" >&2
  exit 1
fi

# Build dirlist.txt with absolute paths
find "$OUT_BASE" -mindepth 1 -maxdepth 2 -type d -print0 \
  | xargs -0 -I{} realpath "{}" > dirlist_PPY.txt

# Safety check: make sure we found something
if [[ ! -s dirlist_PPY.txt ]]; then
  echo "dirlist_PPY.txt is empty — no per-sample output directories found under $OUT_BASE" >&2
  exit 1
fi

echo "Found $(wc -l < dirlist.txt) sample output directories. Saved to dirlist.txt."

# 2) Merge non-reference ME insertions
singularity exec "${sif}" joint_calling_hs \
  -merge_mei \
  -f dirlist_PPY.txt \
  -fa ../data/PPY/assembly/GWHBAOS00000000.genome.fasta.reformat \
  -cohort_name "${SPECIES}" \
  -outdir ../output/PPY/jointcall_out_PlusAdmx  \
  -rep ../data/PPY/TE_annotation_URGI/PPY3_TEannotGr2_refTEs.fa.RMformat

# 3) Merge reference ME polymorphisms
singularity exec /scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif joint_calling_hs \
  -merge_absent_me \
  -f dirlist_PPY.txt \
  -fa ../data/PPY/assembly/GWHBAOS00000000.genome.fasta.reformat \
  -cohort_name "${SPECIES}" \
  -outdir ../output/PPY/jointcall_out_PlusAdmx  \
  -rep ../data/PPY/TE_annotation_URGI/PPY3_TEannotGr2_refTEs.fa.RMformat
  
  
################################################################################
#### (Optional) Step 3. Make a joint call for haplotype phasing

#singularity exec "${sif}" reshape_vcf \
#  -i ../output/PCOM/jointcall_out/"${SPECIES}"_MEI_jointcall.vcf.gz \
#  -a ../output/PCOM/jointcall_out/"${SPECIES}"_MEA_jointcall.vcf.gz \
#  -cohort_name "${SPECIES}"






