#!/usr/bin/env bash

#SBATCH -J s01.run_pixy.sh
#SBATCH -o s01.run_pixy.sh.%J.out
#SBATCH -e s01.run_pixy.sh.%J.err
#SBATCH -c 10
#SBATCH --mem=20G

# @File     :   s01.run_pixy.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/11/03 14:12:40
# Description:
	# 

# Update: 2024-01-17 12:27:20
# 1. define the cores and memory

source s00.config.sh

${LOAD_PIXY}

pixy --stats pi fst dxy \
--vcf "${INPUT_VCF}" \
--populations "${INPUT_POP}" \
--window_size 10000 \
--n_cores 10 \
--output_folder "${OUTPUT_DIR}" \
--output_prefix "${PREFIX}" \
--bypass_invariant_check no

