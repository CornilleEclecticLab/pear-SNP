#!/bin/bash
#SBATCH -J Clone_4.s01.raw_qc
#SBATCH -o Clone_4.s01.raw_qc.out
#SBATCH -e Clone_4.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_4
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Clone_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_2.fq.gz &&\ 
echo done
