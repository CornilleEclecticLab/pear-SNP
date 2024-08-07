#!/usr/bin/env bash

# @File     :   s03.generate_multihetsep
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/11/09 21:14:28
#Description:
	# 


# Set work path
source s00.load_environment.sh


generate_multihetsep.py \
    --mask ${CHR}.${REF_PREFIX}.mask.bed\
	chr76.A301.vcf.gz  chr76.P11-1.vcf.gz\
	> ${WORK_DIR}/output/multihetsep.txt 

