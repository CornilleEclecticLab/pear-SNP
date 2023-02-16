#!/usr/bin/env bash

#SBATCH -J s02.b1.reduce_samples_filter_geno_maf.sh
#SBATCH -o s02.b1.reduce_samples_filter_geno_maf.sh.out
#SBATCH -e s02.b1.reduce_samples_filter_geno_maf.sh.error

# @File     :   s02.b1.reduce_samples_filter_geno_maf.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/02/15 15:44:17
#Description:
	# 

module purge
module load bioinfo/bcftools-1.14

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"



bcftools view $OUTPUT/pear.withLoquat.Combine_Chr.vcf.gz \
    -S ^s02.b2.reduce_samples_filter_geno_maf.LowQualOrCloneRemoveList.txt \
| bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 \
    --threads 8 \
    -o $OUTPUT/pear.withLoquat.Combine_Chr.Reduce3LowQualEuropeWild20Clone.geno20_maf005.vcf.gz 
