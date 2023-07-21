#!/usr/bin/env bash

# @File     :   s00.install_genmap_via_mamba.bash
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/07/18 21:07:35
#Description:
	# 


module load conda


# add bioconda with depended channels and set priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict


# install genmap
WORKDIR="~/work/conda/env/genmap"
mamba create -p $WORKDIR -c bioconda genmap

# init mamba
mamba init bash
source ~/.bashrc