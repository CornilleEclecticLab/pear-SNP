#!/usr/bin/env bash

#SBATCH -J s01.prepear_MEGAnE_k-mer.sh
#SBATCH -o s01.prepear_MEGAnE_k-mer.sh.%J.out
#SBATCH -e s01.prepear_MEGAnE_k-mer.sh.%J.err

# @File     :   s01.prepear_MEGAnE_k-mer.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 17:26:10
#Description:
    # 

source s00.config.sh
module load tmux

singularity exec ${sif} build_kmerset \
-fa genome_file.fasta \
-prefix reference_Pyrifolia_Cuiguan \
-outdir megane_kmer_set
