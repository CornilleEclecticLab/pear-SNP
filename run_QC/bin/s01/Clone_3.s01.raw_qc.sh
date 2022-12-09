#!/bin/bash
#SBATCH -J Clone_3.s01.raw_qc
#SBATCH -o Clone_3.s01.raw_qc.out
#SBATCH -e Clone_3.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_3
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_3 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_3/Clone_3_EKDN220024215-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_3 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_3/Clone_3_EKDN220024215-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_3 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_3/Clone_3_EKDN220024215-1A_H3N2VDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_3 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_3/Clone_3_EKDN220024215-1A_H3N2VDSX5_L1_2.fq.gz &&\ 
echo done
