#!/bin/bash
#SBATCH -J s04.clean_multiqc
#SBATCH -o s04.clean_multiqc.out
#SBATCH -e s04.clean_multiqc.err
module purge

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

module load bioinfo/MultiQC-v1.11

multiqc -n /work/zruilin/Peach/QC/output/s04.clean_multiqc/Peach_clean_report.html \
/work/zruilin/Peach/QC/output/s03.fastp_and_qc.*/*/
