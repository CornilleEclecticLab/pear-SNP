#!/bin/bash
#SBATCH -J P1908.s01.raw_qc
#SBATCH -o P1908.s01.raw_qc.out
#SBATCH -e P1908.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/P1908
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P1908 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H33GNDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P1908 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H33GNDSX5_L2_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P1908 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H3MWGDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/P1908 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H3MWGDSX5_L3_2.fq.gz &&\ 
echo done
