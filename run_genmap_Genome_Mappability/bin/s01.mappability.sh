#!/usr/bin/env bash

#SBATCH -J s01.mappability.sh
#SBATCH -o s01.mappability.NOGIT.%J.out
#SBATCH -e s01.mappability.NOGIT.%J.err
#SBATCH -c 40 


# @File     :   s01.mappability.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/07/18 21:32:38
#Description:
	# 
module purge
module load conda


# activate genmap
source activate ~/work/conda/env/genmap

# genmap map -K 140 -E 0 -T 40 -I ../input/index_genmap -O ../output -t -w -bg 
genmap map -K 140 -E 0 -T 40 -I ../input/index_genmap_Pyrus_communis -O ../outuput/Pyrus_communis -t -w -bg
