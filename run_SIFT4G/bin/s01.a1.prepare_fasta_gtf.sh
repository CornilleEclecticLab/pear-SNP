#!/bin/bash

#SBATCH -J s01.a1.prepare_fasta_gtf.sh
#SBATCH -o s01.a1.prepare_fasta_gtf.sh.%J.out
#SBATCH -e s01.a1.prepare_fasta_gtf.sh.%J.err

# @File     :   s01.a1.prepare_fasta_gtf.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/05/28 18:25:12
#Description:
	# 

module load seqkit
module load gffread
source s00.config.sh


seqkit grep -f ${CHR_LIST} ${GENOME_FASTA_FULL} -o ${GENOME_FASTA_FULL%.fasta}.chr.fa

gzip ${GENOME_FASTA_FULL%.fasta}.chr.fa

grep -f ${CHR_LIST} ${GFF_FILE} > "${GFF_FILE%.gff}.chr.gff"

gffread ${GFF_FILE%.gff}.chr.gff -T -o ${GFF_FILE%.gff}.chr.gtf

gzip ${GFF_FILE%.gff}.chr.gtf
