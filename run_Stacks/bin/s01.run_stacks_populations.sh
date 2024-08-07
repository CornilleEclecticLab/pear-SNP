#!/usr/bin/env bash

#SBATCH -J s01.run_stacks_populations.sh
#SBATCH -o s01.run_stacks_populations.sh.%J.out
#SBATCH -e s01.run_stacks_populations.sh.%J.err

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
    -V ../input/pear_Dec2023.noClone.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf \
    -M ../input/sample_tab_population.txt \
    -O "$OUT_DIR"
