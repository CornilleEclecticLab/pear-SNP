#!/usr/bin/env bash

# @File     :   s03.b1.splite_mask_file_by_chromosome.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/05/22 19:16:52
#Description:
    # 


# Record software version
cd ../input
mkdir mask && cd mask
for i in $(cat ../s02.chromosome_list.txt) ; do echo $i | xargs -I % grep % s03.mask_pass.bed > $i.mask_pass.bed ; done