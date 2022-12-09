#!/bin/bash
#SBATCH -J Clone_5.s01.raw_qc
#SBATCH -o Clone_5.s01.raw_qc.out
#SBATCH -e Clone_5.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_5
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_5/Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_5/Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_2.fq.gz &&\ 
echo done
