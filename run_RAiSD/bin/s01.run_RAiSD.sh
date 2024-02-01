#!/usr/bin/env bash

# @File     :   s01.run_RAiSD.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/11/17 21:44:21
#Description:
    # 

source s00.load_RAiSD.sh

RAiSD \
    -n EuropeWild \
    -I pear_Jul2023.Combine_Chr.EW.chr76.maf001.vcf \
    -w 50 \
    -M 3 \
    -y 2 \
    -R \
    -f \
    -k 0.05\
    -P

RAiSD \
    -n EuropeWild_Grid \
    -I pear_Jul2023.Combine_Chr.EW.chr76.maf001.vcf \
    -G 9085\
    -R \
    -f \
    -k 0.05\
    -P