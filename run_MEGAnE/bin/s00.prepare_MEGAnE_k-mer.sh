#!/usr/bin/env bash

#SBATCH -J s00.prepare_MEGAnE_k-mer.sh
#SBATCH -o s00.prepare_MEGAnE_k-mer.sh.%J.out
#SBATCH -e s00.prepare_MEGAnE_k-mer.sh.%J.err
#SBATCH -c 1
#SBATCH --mem=20G

# @File     :   s00.prepear_MEGAnE_k-mer.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 17:26:10
#Description:
    # 

source s00.config.sh
module load tmux


singularity exec ${sif} build_kmerset \
-fa ../input/ref_genome.fa \
-prefix reference_Pyrifolia_Cuiguan \
-outdir ../input/megane_kmer_set
