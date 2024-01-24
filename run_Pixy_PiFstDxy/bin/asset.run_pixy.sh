#!/usr/bin/env bash

# @File     :   building.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/05/04 22:53:00
# Description:
	# 
# Update: 2024-01-17 12:27:20
# ChangeLog:
    # 1. using source activate to load pixy conda environment


source activate pixy

pixy --stats pi fst dxy \
--vcf test.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 4 \
--output_folder output \
--output_prefix pixy_output \
--bypass_invariant_check no # Bypass the check for invariant sites. Use with caution!


