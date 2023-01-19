#!/usr/bin/env bash
#SBATCH -J index_samtools
#SBATCH -o index_samtools.out
#SBATCH -e index_samtools.err


module purge
module load bioinfo/samtools-1.14

samtools faidx ../input/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.fasta
