#!/usr/bin/env bash
#SBATCH -J index_samtools
#SBATCH -o index_samtools.out
#SBATCH -e index_samtools.err


module purge

module load samtools/1.14

samtools faidx \
    /shared/ifbstor1/projects/pear_snp3/pear/run_GATK_variant_calling_refer_on_communis/input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta