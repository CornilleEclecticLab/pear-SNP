#!/usr/bin/env bash

#SBATCH -J s03.joint_calling.sh
#SBATCH -o s03.joint_calling.sh.%J.out
#SBATCH -e s03.joint_calling.sh.%J.err

# @File     :   s03.joint_calling.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 18:10:40
#Description:
    # 

source

# first, list up samples (output directories from Step 02) you are going to merge
ls -d /path/to/[all_output_directories] > dirlist.txt

# merge non-reference ME insertions
singularity exec ${sif} joint_calling_hs \
-merge_mei \
-f dirlist.txt \
-fa genome_file.fasta \
-cohort_name test

# merge reference ME polymorphisms
singularity exec ${sif} joint_calling_hs \
-merge_absent_me \
-f dirlist.txt \
-fa genome_file.fasta \
-cohort_name test