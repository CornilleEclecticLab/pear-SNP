#!/bin/bash

#SBATCH -J s02.run_blastp.sh
#SBATCH -o s02.run_blastp.sh.%J.out
#SBATCH -e s02.run_blastp.sh.%J.err
#SBATCH -c 24

# @File     :   s02.run_blastp.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/02/15 17:08:29
#Description:
    # 

source s00.load_blast.sh
source s00_config.py

mkdir -p ${s02_dir}

# blastp
for query in ${s01_dir}/*.fa ; do
    for db in ${db_dir}/*.pdb ; do
        query_name=$(basename ${query} .fa)
        db_name=$(basename ${db} .pdb)
        echo "Running blastp for ${query_name} against ${db_name} ..."
        blastp \
            -query ${query} \
            -db ${db_dir}/${db_name} \
            -outfmt 6 \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads 4 \
            -out ${s02_dir}/${query_name}.${db_name}.blastp.tsv &
    done
done

wait
