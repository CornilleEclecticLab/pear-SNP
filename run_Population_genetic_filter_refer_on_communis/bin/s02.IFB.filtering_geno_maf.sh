#!/bin/bash
# Adaptation by Yuqi NIE on 2022/11/02
# Edited by Yuqi NIE on 2023/08/02
# Updated by Yuqi NIE on 2024-01-05
# Updated by Yuqi NIE on 2024-01-09

#SBATCH -J s02.filtering_geno_maf
#SBATCH -o s02.IFB.remove_clone_filtering_geno_maf.%J.out
#SBATCH -e s02.IFB.remove_clone_filtering_geno_maf.%J.err
#SBATCH -c 8
module purge
module load bcftools/1.14

source s00.config.sh

mkdir -p $OUTPUT

# -S ^list
bcftools view \
-S s02.select_list.txt \
    "${IN_OUTPUT}/${IN_PREFIX}.Combine_chr.vcf.gz" \
    --threads 8 \
| bcftools filter \
    -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' \
    -O z4 \
    --threads 8 \
    -o $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.vcf.gz 
