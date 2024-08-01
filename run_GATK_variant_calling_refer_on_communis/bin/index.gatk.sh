#!/bin/bash
#SBATCH -J index_gatk
#SBATCH -o index_gatk.out
#SBATCH -e index_gatk.err

module purge

module load gatk4/4.1.9.0

gatk CreateSequenceDictionary \
    -R /shared/ifbstor1/projects/pear_snp3/pear/run_GATK_variant_calling_refer_on_communis/input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta
    