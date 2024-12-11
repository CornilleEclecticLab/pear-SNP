#!/bin/bash
#SBATCH --job-name=mkt_analysis
#SBATCH --output=s02.run_MKT_slurm.%J.out 
#SBATCH --error=s02.run_MKT_slurm.%J.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Load R, and parameters
source s00.config.sh

for N in {1..17}; do
  chr="Chr${N}"
  echo $chr
  vcf_file="${vcf_prefix}${chr}${vcf_suffix}"
  Rscript s02.x1.perform_MKT_with_args.R ${vcf_file} ${gff_file} ${fasta_file} ${sample_file} ${outgroup_pop}
done
