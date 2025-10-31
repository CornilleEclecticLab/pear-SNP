#!/usr/bin/env bash

#SBATCH -J s01.run_cd-hit.sh
#SBATCH -o s01.run_cd-hit.sh.%J.out
#SBATCH -e s01.run_cd-hit.sh.%J.err
#SBATCH -c 36
#SBATCH --mem=24G

# @File     :   s01.run_cd-hit.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/03 15:42:52
#Description:
    # 

module load cd-hit/4.8.1

cd-hit-est \
    -i ../input/Inpactor2_library.fasta \
    -T 0 \
    -M 0 \
    -o ../output/reduced_Inpactor2_library_c10.fasta \
    -d 0 \
    -aS 0.8 \
    -c 0.8 \
    -G 0 \
    -g 1 \
    -b 500
