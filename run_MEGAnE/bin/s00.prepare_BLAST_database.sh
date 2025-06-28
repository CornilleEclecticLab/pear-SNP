#!/usr/bin/env bash

#SBATCH -J s00.prepare_BLAST_database.sh
#SBATCH -o s00.prepare_BLAST_database.sh.%J.out
#SBATCH -e s00.prepare_BLAST_database.sh.%J.err

# @File     :   s00.prepare_BLAST_database.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/29 13:13:22
#Description:
    # 

source s00_config.py

mkdir -p ../input/ref_genome_blastdb

singularity exec $sif /usr/local/bin/ncbi-blast-2.12.0+/bin/makeblastdb \
-in ../input/${ref_genome_fa} \
-out ../input/ref_genome_blastdb/${ref_genome_fa} \
-dbtype nucl \
-parse_seqids