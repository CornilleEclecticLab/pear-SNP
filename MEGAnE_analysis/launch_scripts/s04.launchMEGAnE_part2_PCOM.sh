#!/bin/bash
#SBATCH --job-name=PCOM_MEGAnE_part2_admixed
#SBATCH --output=../logs/MEGAnE_part2_PCOM_admixed.out
#SBATCH --error=../logs/MEGAnE_part2_PCOM_admixed.err
#SBATCH --time=6-00:00:00
#SBATCH -p compute
#SBATCH -c 64
#SBATCH --mem=150G

set -euo pipefail
# Environment setup
module load singularity
module load all gencore/3
module load blast/2.16.0
module load repeatmasker/4.1.9

# Singularity image for MEGAnE
sif=/scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif

################################################################################
# Step 1. Call and genotype polymorphic MEs
# Required flags for call_genotype (from MEGAnE docs):
#  -repout str      RepeatMasker output masked using the input RepBase file
#  -repremove str   File listing non-ME repeat classes
#  -pA_ME str       File listing ME classes with polyA tail
#  -mainchr str     File listing main chromosomes for non-human genomes

# Root folder containing the CRAM files (no nested structure assumptions)
ROOT="../data/PCOM/CRAM" #for admixed the source folder was /scratch/yn2515/pear_snp_tip/cram_files_with_Admix/western_pears_ref_PCOM

# Find all CRAM files under ROOT
find "$ROOT" -type f -name '*.cram' -print0 | while IFS= read -r -d '' cram; do
  file_name="$(basename "$cram")"

  # Extract SAMPLE as the 3rd field when splitting by '.' (e.g., pear.Teng.P6-4.sorted.markdu.cram -> P6-4)
  # If the filename doesn't have at least 3 fields, skip safely.
  IFS='.' read -r f1 f2 SAMPLE rest <<< "$file_name" || true
  if [[ -z "${SAMPLE:-}" ]]; then
    echo "⚠️  Skipping (cannot parse SAMPLE as 3rd dot-separated field): $cram" >&2
    continue
  fi

  # Output directory per-sample
  # outdir="../output/PCOM/MEGAnE_result/$SAMPLE"
  outdir="../output/PCOM/MEGAnE_result_admixed/$SAMPLE"

  # Skip if already processed (directory exists and not empty)
  if [[ -d "$outdir" && -n "$(ls -A "$outdir" 2>/dev/null)" ]]; then
    echo "⏩  Already processed: SAMPLE='$SAMPLE' (outdir exists and is not empty)"
    continue
  fi

  mkdir -p "$outdir"

  echo "▶️  Processing: SAMPLE='$SAMPLE' | CRAM='$cram' | OUTDIR='$outdir'"

  # Run MEGAnE call_genotype inside the Singularity container
  singularity exec "${sif}" call_genotype \
    -i "$cram" \
    -fa ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat \
    -fadb ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat \
    -mk ../output/PCOM/megane_kmer_set/reference_PCOM_genome.mk \
    -outdir "$outdir" \
    -p 128 \
    -rep ../data/PCOM/TE_annotation_URGI/PCOM_TEannotGr2_refTEs.fa.RMformat \
    -repout ../data/PCOM/TE_annotation_URGI/PCOM_TE_annotation.gff3.out.reformat \
    -repremove ../data/PCOM/non_ME_rep.txt \
    -pA_ME ../data/PCOM/ME_with_pA.txt \
    -mainchr ../data/PCOM/main_chr_list.txt \
    -no_sex_chr
done
