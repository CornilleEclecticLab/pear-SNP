#!/bin/bash
#SBATCH -J s02.raw_multiqc
#SBATCH -o s02.raw_multiqc.out
#SBATCH -e s02.raw_multiqc.err
module purge

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

module load bioinfo/MultiQC-v1.11

multiqc -n /work/zruilin/Peach/QC/output/s02.raw_multiqc/pear_raw_report.html \
/work/zruilin/Peach/QC/output/s01.raw_qc/*/

