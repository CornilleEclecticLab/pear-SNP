#!/bin/bash

# @File     :   s05.a1.run_regioneR_permutation_test.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/07/28 09:54:02
#Description:
	# 



# Count *.genes.gff files



# Export the list
ls ../output/s04.catch_outliers/Pixy_TajimaD.0_0*.gff > s05.a1.gff_list.txt

N=$(wc -l < s05.a1.gff_list.txt)

# Submit array
sbatch --array=1-$N%200 s05.x1.regioneR_array.sh