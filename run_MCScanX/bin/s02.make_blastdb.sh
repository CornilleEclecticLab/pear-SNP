#!/usr/bin/env bash

#SBATCH -J s02.make_blastdb.sh
#SBATCH -o s02.make_blastdb.sh.%J.out
#SBATCH -e s02.make_blastdb.sh.%J.err

# @File     :   s02.make_blastdb.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/02/10 12:00:42
#Description:
    # 

source s00.load_blast.sh
source s00_config.py

mkdir -p $s02_dir

for i in ${s01_dir}/*.fa; do
    echo "Processing file: $i"
    base=$(basename "$i" ".fa")
    echo "Base name: $base"

    makeblastdb -in "$i" \
        -dbtype prot \
        -parse_seqids \
        -out "${s02_dir}/${base}"
done