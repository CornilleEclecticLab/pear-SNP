#!/bin/bash
#SBATCH --job-name=PCOM_MEGAnE_part3
#SBATCH --output=../logs/PCOM_MEGAnE_output_part3.txt
#SBATCH --error=../logs/PCOM_MEGAnE_error_part3.txt
#SBATCH --time=6-00:00:00
#SBATCH -p compute
#SBATCH -c 32
#SBATCH --mem=64G

set -euo pipefail

module load singularity

sif=/scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif

# Base directory where per-sample outputs were written in Step 1
OUT_BASE="../output/PCOM/MEGAnE_result"

# Species/cohort name variable
SPECIES='PCOM'

################################################################################
#### Step 2. Joint calling
# first, list up samples (output directories from Step 1) you are going to merge
#ls -d /path/to/[all_output_directories] > dirlist.txt

# 1) List all per-sample MEGAnE output directories (one level under OUT_BASE)
#    Write ABSOLUTE paths into dirlist.txt (some tools prefer absolute paths).
#    If you want to be stricter, you can filter only directories containing expected result files.
if [[ ! -d "$OUT_BASE" ]]; then
  echo "No output base directory found at $OUT_BASE" >&2
  exit 1
fi

# Build dirlist.txt with absolute paths
find "$OUT_BASE" -mindepth 1 -maxdepth 2 -type d -print0 \
  | xargs -0 -I{} realpath "{}" > dirlist.txt

# Safety check: make sure we found something
if [[ ! -s dirlist.txt ]]; then
  echo "dirlist.txt is empty — no per-sample output directories found under $OUT_BASE" >&2
  exit 1
fi

echo "Found $(wc -l < dirlist.txt) sample output directories. Saved to dirlist.txt."

# 2) Merge non-reference ME insertions
singularity exec "${sif}" joint_calling_hs \
  -merge_mei \
  -f dirlist.txt \
  -fa ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat \
  -cohort_name "${SPECIES}" \
  -outdir ../output/PCOM/jointcall_out_PlusAdmx  \
  -rep ../data/PCOM/TE_annotation_URGI/PCOM_TEannotGr2_refTEs.fa.RMformat

# 3) Merge reference ME polymorphisms
singularity exec /scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif joint_calling_hs \
  -merge_absent_me \
  -f dirlist.txt \
  -fa ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat \
  -cohort_name "${SPECIES}" \
  -outdir ../output/PCOM/jointcall_out_PlusAdmx  \
  -rep ../data/PCOM/TE_annotation_URGI/PCOM_TEannotGr2_refTEs.fa.RMformat
  
  
################################################################################
#### (Optional) Step 3. Make a joint call for haplotype phasing

#singularity exec "${sif}" reshape_vcf \
#  -i ../output/PCOM/jointcall_out/"${SPECIES}"_MEI_jointcall.vcf.gz \
#  -a ../output/PCOM/jointcall_out/"${SPECIES}"_MEA_jointcall.vcf.gz \
#  -cohort_name "${SPECIES}"






