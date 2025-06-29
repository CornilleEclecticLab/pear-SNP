#!/bin/bash

#SBATCH -J s02.x2.reshape_vcf.sh
#SBATCH -o s02.x2.reshape_vcf.sh.log/%J.out
#SBATCH -e s02.x2.reshape_vcf.sh.log/%J.err

# @File     :   s03.x2.reshape_vcf.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/06/19 14:53:18
#Description:
	# 

source s00_config.py
source s00.load_bcftools.sh

cohort_name=$1

temp_dir="../tmp"
mkdir -p "${temp_dir}"

singularity exec ${sif} reshape_vcf \
	-i "${s02_output_dir}/${cohort_name}/${cohort_name}_MEI_jointcall.vcf.gz" \
	-a "${s02_output_dir}/${cohort_name}/${cohort_name}_MEA_jointcall.vcf.gz" \
	-cohort_name "${cohort_name}" \
	-outdir "${s02_output_dir}/${cohort_name}"

bcftools view \
	"${s02_output_dir}/${cohort_name}/${cohort_name}_biallelic.vcf.gz" \
	--apply-filter "PASS" \
| bcftools sort \
    --temp-dir "${temp_dir}" \
	--output-type z \
	--output-file "${s02_output_dir}/${cohort_name}/${cohort_name}_biallelic_pass.vcf.gz"

tabix -p vcf \
	"${s02_output_dir}/${cohort_name}/${cohort_name}_biallelic_pass.vcf.gz"