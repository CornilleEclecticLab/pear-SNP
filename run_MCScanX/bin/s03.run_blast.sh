#!/usr/bin/env bash

#SBATCH -J s03.run_blast.sh
#SBATCH -o s03.run_blast.sh.%J.out
#SBATCH -e s03.run_blast.sh.%J.err
#SBATCH -c 40
#SBATCH --array=0-3

# @File     :   s03.run_blast.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/02/10 12:23:26
#Description:
    # 

source s00.load_blast.sh
source s00_config.py


comm_db=$(basename "${comm_pep_fa}" ".fasta").formatted
pyri_db=$(basename "${pyri_pep_fa}" ".fasta").formatted

comm_pep_formatted_fa="${comm_db}.fa"
pyri_pep_formatted_fa="${pyri_db}.fa"

# Define arrays
query=( "${comm_pep_formatted_fa}" "${pyri_pep_formatted_fa}" "${comm_pep_formatted_fa}" "${pyri_pep_formatted_fa}" )
db=( "${comm_db}" "${comm_db}" "${pyri_db}" "${pyri_db}" )
out=( "comm_comm" "pyri_comm" "comm_pyri" "pyri_pyri" )

index=${SLURM_ARRAY_TASK_ID}

echo "Running blastp on ${query[$index]} against ${db[$index]}"

blastp -query "${s01_dir}/${query[$index]}" \
    -db "${s02_dir}/${db[$index]}" \
    -out "${s03_dir}/${out[$index]}.blast.txt" \
    -outfmt 6 \
    -num_threads 40 \
    -evalue 1e-10 \
    -max_target_seqs 5
