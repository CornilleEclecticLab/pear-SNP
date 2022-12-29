#!/bin/bash
#SBATCH -J s04.clean_multiqc
#SBATCH -o s04.clean_multiqc.out
#SBATCH -e s04.clean_multiqc.err
module purge

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

multiqc -n /data/atipe-workspace/ynie/pear/run_QC/output/s04.clean_multiqc/Pear_clean_report.html \
/data/atipe-workspace/ynie/pear/run_QC/output/s03.fastp_and_qc/*/*/
