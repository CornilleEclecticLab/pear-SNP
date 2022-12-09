#!/bin/bash
#SBATCH -J S7190.s01.raw_qc
#SBATCH -o S7190.s01.raw_qc.out
#SBATCH -e S7190.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/S7190
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S7190 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S7190 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S7190 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H3MWGDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S7190 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H3MWGDSX5_L3_2.fq.gz &&\ 
echo done
