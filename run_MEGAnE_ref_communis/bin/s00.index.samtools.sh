#!/usr/bin/env bash
#SBATCH -J index_samtools
#SBATCH -o s00.index_samtools.%J.out
#SBATCH -e s00.index_samtools.%j.err


module purge
module load samtools/1.14

samtools faidx ../input/ref_genome.fa
