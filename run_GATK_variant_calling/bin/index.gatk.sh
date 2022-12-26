#!/bin/bash
#SBATCH -J index_gatk
#SBATCH -o index_gatk.out
#SBATCH -e index_gatk.err

module purge

module load bioinfo/gatk-4.1.7.0

gatk CreateSequenceDictionary -R ../input/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.fasta -O ../input/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.dict

