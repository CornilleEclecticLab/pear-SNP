#!/bin/bash

#SBATCH --job-name=s05.regioneR
#SBATCH --output=s05.x1_log/s05.run_regioneR.%J.out
#SBATCH --error=s05.x1_log/s05.run_regioneR.%J.err
#SBATCH -c 2
#SBATCH --time=00:20:00



# @File     :   s05.run_regioneR.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/07/27 14:53:47
#Description:
	# 

module load r/4.3.1


GFF_GENOME="../input/GWHBAOS00000000.chr_list.gene.slop1kb.gff"

GENOME_DEF="../input/karyotype.txt"

GFF_LIST="s05.a1.gff_list.txt"

OUTPUT_DIR="../output/s05.regioneR_permutation_test"

mkdir -p ${OUTPUT_DIR}

GFF_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $GFF_LIST)

echo "Processing GFF file: ${GFF_FILE}"
Rscript s05.x2.regioneR_permutation_test.R \
	${GFF_FILE} \
	${GFF_GENOME} \
	${GENOME_DEF} \
	${OUTPUT_DIR}