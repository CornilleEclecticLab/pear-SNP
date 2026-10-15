#!/bin/bash

#SBATCH -J s08.b2.enrichment_pi_matched.slurm.sh
#SBATCH -o s08.b2.enrichment_pi_matched.slurm.sh.%J.out
#SBATCH -e s08.b2.enrichment_pi_matched.slurm.sh.%J.err
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=04:00:00


module load r/4.3.1
export TMPDIR=$(realpath ../tempdir)   # do not use TMPDIR from the submitting session, as it was full

Rscript s08.x1.enrichment_pi_matched.R \
  --lineage Oriental \
  --dinv "../output/Dinvestigate/*/*_localFstats__50_25.txt" \
  --sweep_bed "../input/selection_bed/selection.merge_across_methods.*.merged.bed" \
  --mask_bed ../input/GWHBAOS00000000.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed \
  --genes ../input/GWHBAOS00000000.chr_list.gene.slop1kb.gff \
  --pixy_pi "../../run_Pixy_PiFstDxy/output/s02.run_pixy_by_chromosome_for_local_D/*_pi.txt" \
  --pi_mode bed \
  --s06_summary_dir \
  ../output/Dinvestigate_summary \
  --threads 16
