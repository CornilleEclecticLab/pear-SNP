#!/bin/bash
#SBATCH -J S5398.s01.raw_qc
#SBATCH -o S5398.s01.raw_qc.out
#SBATCH -e S5398.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/S5398
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5398 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5398/S5398_EKDN220024243-1A_H3MGWDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5398 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5398/S5398_EKDN220024243-1A_H3MGWDSX5_L3_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5398 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5398/S5398_EKDN220024243-1A_H3TGLDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5398 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5398/S5398_EKDN220024243-1A_H3TGLDSX5_L3_2.fq.gz &&\ 
echo done
