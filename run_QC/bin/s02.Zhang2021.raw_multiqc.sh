#!/bin/bash
#SBATCH -J s02.raw_multiqc
#SBATCH -o s02.raw_multiqc.out
#SBATCH -e s02.raw_multiqc.err

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8


multiqc -n ../output/s02.raw_multiqc/pear_Zhang2021_report.html \
../output/s01.raw_qc/Pear_zhang2021/*/

