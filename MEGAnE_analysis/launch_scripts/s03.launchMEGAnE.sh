#!/bin/bash
#SBATCH --job-name=MEGAnE_part1
#SBATCH --output=../logs/MEGAnE_part1.out
#SBATCH --error=../logs/MEGAnE_part1.err
#SBATCH --time=3:00:00
#SBATCH -p compute
#SBATCH -c 8
#SBATCH --mem=16G

set -euo pipefail

module purge
module load singularity
module load all gencore/3
module load blast/2.16.0
module load repeatmasker/4.1.9
module load samtools

sif=/scratch/ss20440/Aphid_project/MEGAnE/MEGAnE.sif

mkdir -p ../logs
mkdir -p ../output

################################################################################
#### TORM Reformat fasta genome files, tig00001236_quiver_pilon_obj becomes tig00001236
# TORM awk '/^>/{sub(/_quiver.*/,"")}1' \
# TORM ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta \
# TORM  > ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat

################################################################################
#### Step 0. Prepare MEGAnE k-mer file
singularity exec ${sif} build_kmerset \
 -fa ../data/PPY/assembly/GWHBAOS00000000.genome.fasta.reformat \
 -prefix reference_PPY_genome \
 -outdir ../output/PPY/megane_kmer_set

# singularity exec ${sif} build_kmerset \
# -fa ../data/TE_annotation_URGI/PCOM/PyrusCommunis_BartlettDHv2.0.fasta.reformat \
# -prefix reference_PCOM_genome \
# -outdir ../output/PCOM/megane_kmer_set

################################################################################
#### Step 1.0 Prepare fadb file
makeblastdb -in ../data/PPY/assembly/GWHBAOS00000000.genome.fasta.reformat -dbtype nucl
# makeblastdb -in ../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat -dbtype nucl


