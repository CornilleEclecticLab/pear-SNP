#!/bin/bash

#SBATCH -J s01.run_orthofinder.sh
#SBATCH -o s01.run_orthofinder.sh.%J.out
#SBATCH -e s01.run_orthofinder.sh.%J.err
#SBATCH -c 32

# @File     :   s01.run_orthofinder.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/05/05 10:46:04
#Description:
    # 

source s00.load_orthofinder.sh

orthofinder -f ../input \
            -o ../output 
