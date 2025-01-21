#!/usr/bin/env bash

#SBATCH -J s01.run_count_genes_by_windows_from_gff.sh
#SBATCH -o s01.run_count_genes_by_windows_from_gff.sh.%J.out
#SBATCH -e s01.run_count_genes_by_windows_from_gff.sh.%J.err
#SBATCH -c 1
#SBATCH --mem=10G

# @File     :   s01.run_count_genes_by_windows_from_gff.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/01/15 10:49:32
#Description:
    # 

module load python

python3 s01.count_genes_by_windows_from_gff.py
