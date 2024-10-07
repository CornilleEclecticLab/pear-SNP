#!/usr/bin/env bash

#SBATCH -J s01.run_stacks_populations.sh
#SBATCH -o s01.run_stacks_populations.sh.%J.out
#SBATCH -e s01.run_stacks_populations.sh.%J.err
#SBATCH -c 8
#SBATCH --mem=64G

# @File     :   s01.run_stacks_populations.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/02/23 17:46:07
#Description:
    # 

module load gcc/11.2.0
source activate stacks

OUT_DIR="../output/s01.run_stacks_populations"
mkdir -p  "${OUT_DIR}"

populations \
    -k \
    -t 8 \
    -V ../input/input.vcf.gz \
    -M ../input/sample_tab_population.txt \
    -O "$OUT_DIR"
