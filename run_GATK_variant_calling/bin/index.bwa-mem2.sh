#!/bin/bash
#SBATCH -J index_bwa-mem2
#SBATCH -o index_bwa-mem2.out
#SBATCH -e index_bwa-mem2.err

module purge
module load bioinfo/bwa-mem2-2.2

bwa-mem2.sse41 index ../input/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.fasta
