#!/usr/bin/env bash
#SBATCH -J index_samtools
#SBATCH -o index_samtools.out
#SBATCH -e index_samtools.err


module purge
module load samtools/1.14

samtools faidx ../input/ref_genome.fa
