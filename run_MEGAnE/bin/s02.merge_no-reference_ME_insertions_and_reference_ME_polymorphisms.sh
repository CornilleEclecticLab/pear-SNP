#!/usr/bin/env bash

#SBATCH -J s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh
#SBATCH -o s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh.%J.out
#SBATCH -e s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh.%J.err

# @File     :   s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/05 21:00:38
#Description:
    # 



# Usage:
# List up samples (output directories from Step 1) of each cohort, 
# And then merge non-reference ME insertions and reference ME polymorphisms for each cohort.
: <<'EOF'
for path in ../output/s01.*/*; do
    popbase=$(basename "$path")
    realpath "$path"/* > "../input/s02.pop_path.$popbase.txt"
    sbatch s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh "../input/s02.pop_path.$popbase.txt" "$popbase" 
done
EOF

source s00.config.sh


input_dir="../input"
output_dir="../output/s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms"
mkdir -p "${output_dir}"

# Merge non-reference ME insertions
singularity exec "${sif}" joint_calling_hs \
    -merge_mei \
    -f "$1" \
    -fa "${input_dir}/ref_genome.fa" \
    -rep "${input_dir}/reduced_Inpactor2_library_c10.fasta" \
    -cohort_name "$2" \
    -outdir "${output_dir}"

# Merge reference ME polymorphisms
singularity exec "${sif}" joint_calling_hs \
    -merge_absent_me \
    -f "$1" \
    -fa "${input_dir}/ref_genome.fa" \
    -rep "${input_dir}/reduced_Inpactor2_library_c10.fasta" \
    -cohort_name "$2" \
    -outdir "${output_dir}"