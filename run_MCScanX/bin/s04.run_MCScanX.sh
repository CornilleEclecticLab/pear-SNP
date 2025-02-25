#!/usr/bin/env bash

#SBATCH -J s04.run_MCScanX.sh
#SBATCH -o s04.run_MCScanX.sh.%J.out
#SBATCH -e s04.run_MCScanX.sh.%J.err

# @File     :   s04.run_MCScanX.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/02/10 15:14:11
#Description:
    # 

source s00.load_MCScanX.sh
source s00_config.py

mkdir -p $s04_dir

cat ${s03_dir}/pyri_comm.blast.txt > ${s04_dir}/pyri_comm.blast
cat ${s03_dir}/comm_pyri.blast.txt > ${s04_dir}/comm_pyri.blast

MCScanX ${s04_dir}/pyri_comm
MCScanX ${s04_dir}/comm_pyri