#!/bin/bash
#SBATCH -J S5543.s01.raw_qc
#SBATCH -o S5543.s01.raw_qc.out
#SBATCH -e S5543.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/S5543
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5543 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5543/S5543_EKDN220024246-1A_H3MGWDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5543 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5543/S5543_EKDN220024246-1A_H3MGWDSX5_L3_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5543 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5543/S5543_EKDN220024246-1A_H3NGMDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S5543 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5543/S5543_EKDN220024246-1A_H3NGMDSX5_L3_2.fq.gz &&\ 
echo done
