#!/usr/bin/env bash

# @File     :   s00.download_tontainer_for_smcpp.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/07/05 13:47:48
#Description:
	# 


# load singularity modules
module load system/singularity-3.7.3

# download container
singularity pull smcpp.sif docker://terhorst/smcpp:latest

singularity run smcpp.sif version > s00.smcpp_version.txt
