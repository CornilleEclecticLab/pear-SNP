#!/bin/bash
#SBATCH -J index_bwa-mem2
#SBATCH -o index_bwa-mem2.out
#SBATCH -e index_bwa-mem2.err

module load conda

source /shared/ifbstor1/home/ynie/.bashrc

conda activate bwa-mem2-2.2

bwa-mem2.sse41 index \
    /shared/ifbstor1/projects/pear_snp3/pear/run_GATK_variant_calling_refer_on_communis/input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta
