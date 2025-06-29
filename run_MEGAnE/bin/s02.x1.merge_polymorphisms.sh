#!/bin/bash

#SBATCH -J s02.x1.joint_calling.log.sh
#SBATCH -o s02.x1.merge_polymorphisms_joint_calling.log/%J.out
#SBATCH -e s02.x1.merge_polymorphisms_joint_calling.log/%J.err

# @File     :   s02.x1.merge_polymorphisms.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/05 21:00:38
#Description:
    #
    # This script is expected to be called by the main script:
    # s02.launch_joint_calling.sh

source s00_config.py

input_dir="../input"
output_dir="${s02_output_dir}/$2"
mkdir -p "${output_dir}"

# Merge non-reference ME insertions
singularity exec "${sif}" joint_calling_hs \
    -merge_mei \
    -f "$1" \
    -fa "${input_dir}/${ref_genome_fa}" \
    -rep "${input_dir}/${rep_library_fa}" \
    -cohort_name "$2" \
    -outdir "${output_dir}"

# Merge reference ME polymorphisms
singularity exec "${sif}" joint_calling_hs \
    -merge_absent_me \
    -f "$1" \
    -fa "${input_dir}/${ref_genome_fa}" \
    -rep "${input_dir}/${rep_library_fa}" \
    -cohort_name "$2" \
    -outdir "${output_dir}"