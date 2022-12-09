#!/bin/bash
#SBATCH -J S8315.s01.raw_qc
#SBATCH -o S8315.s01.raw_qc.out
#SBATCH -e S8315.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/S8315
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S8315 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8315/S8315_EKDN220024259-1A_H3MGWDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S8315 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8315/S8315_EKDN220024259-1A_H3MGWDSX5_L2_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S8315 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8315/S8315_EKDN220024259-1A_H3TGLDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S8315 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8315/S8315_EKDN220024259-1A_H3TGLDSX5_L3_2.fq.gz &&\ 
echo done
