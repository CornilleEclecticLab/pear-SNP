#!/usr/bin/env bash

#SBATCH -J s01.call_and_genotype_polymorphic_MEs.sh
#SBATCH -o s01.call_and_genotype_polymorphic_MEs.sh.%J.out
#SBATCH -e s01.call_and_genotype_polymorphic_MEs.sh.%J.err
#SBATCH -c 16          # Booked 16 but used 4 cores during testing. Job ID: 38787863
#SBATCH --mem=120G     # Booked 120G but used 39.36G during testing. Time 04:18:04

# @File     :   s01.call_and_genotype_polymorphic_MEs.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 18:03:22
#Description:
    # 

source s00.config.sh
input_dir="../input"

outdir='../output/CPP1_No_Skip_Unmapped'
mkdir -p ${outdir}

# Run with a custom repeat library and annotation
singularity exec ${sif} call_genotype \
    -i ${input_dir}/cram_files/pear.CPP_JPP_WPU.CPP1.sorted.markdu.cram \
    -fa ${input_dir}/ref_genome.fa \
    -mk ${input_dir}/megane_kmer_set/reference_Pyrifolia_Cuiguan.mk \
    -fadb ${input_dir}/ref_genome_blastdb/ref_genome_blastdb \
    -rep ${input_dir}/reduced_Inpactor2_library_c10.fasta \
    -repout ${input_dir}/GWHBAOS00000000.genome.fasta.out \
    -repremove ${input_dir}/non_ME_rep.blank.txt \
    -pA_ME ${input_dir}/ME_with_pA.blank.txt \
    -mainchr ${input_dir}/main_chr_list.txt \
    -sample_name CPP1_NSU \
    -no_sex_chr \
    -outdir ${outdir} \
    -p 16 
    # -skip_unmapped 


# Two tips:
# 1. The depth is important,
# 2. The reads should use the same sequencing platform.
