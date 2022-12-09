#!/bin/bash
#SBATCH -J ch154_1.s01.raw_qc
#SBATCH -o ch154_1.s01.raw_qc.out
#SBATCH -e ch154_1.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/ch154_1
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch154_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H33GNDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch154_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H33GNDSX5_L2_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch154_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch154_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_2.fq.gz &&\ 
echo done
