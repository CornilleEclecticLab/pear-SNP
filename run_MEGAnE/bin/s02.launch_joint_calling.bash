#!/bin/bash

# @File     :   s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh
# @Version  :   2.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/05 21:00:38
#Description:
    # 


# Usage:
# bash s02.joint_calling.sh



# List up samples (output directories from Step 1) of each cohort, 
# And then merge non-reference ME insertions and reference ME polymorphisms for each cohort.

: > "../input/s02.pop_path.PPY_joint.txt"

for path in ../output/s01.*/*; do
    popbase=$(basename "$path")
    realpath "$path"/* > "../input/s02.pop_path.$popbase.txt"
    realpath "$path"/* >> "../input/s02.pop_path.PPY_joint.txt"
    
    # Joint calling for each population
    # slurm_jobID=$(sbatch --parsable s02.x1.merge_polymorphisms.sh "../input/s02.pop_path.$popbase.txt" "$popbase")
    # sbatch --dependency=afterok:$slurm_jobID s02.x2.reshape_vcf.sh "$popbase"
    # sbatch s02.x2.reshape_vcf.sh "$popbase"
done

# Job for joint all populations
slurm_jobID=$(sbatch --parsable  s02.x1.merge_polymorphisms.sh  ../input/s02.pop_path.PPY_joint.txt  PPY_joint)
sbatch --dependency=afterok:$slurm_jobID s02.x2.reshape_vcf.sh  PPY_joint

