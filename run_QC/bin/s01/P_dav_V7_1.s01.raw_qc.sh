#!/bin/bash
#SBATCH -J P_dav_V7_1.s01.raw_qc
#SBATCH -o P_dav_V7_1.s01.raw_qc.out
#SBATCH -e P_dav_V7_1.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/P_dav_V7_1
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P_dav_V7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P_dav_V7_1/P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P_dav_V7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P_dav_V7_1/P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_2.fq.gz &&\ 
echo done
