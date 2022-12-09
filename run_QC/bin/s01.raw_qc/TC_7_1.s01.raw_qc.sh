#!/bin/bash
#SBATCH -J TC_7_1.s01.raw_qc
#SBATCH -o TC_7_1.s01.raw_qc.out
#SBATCH -e TC_7_1.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/TC_7_1
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_7_1/TC_7_1_EKDN220024264-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_7_1/TC_7_1_EKDN220024264-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_7_1/TC_7_1_EKDN220024264-1A_H3MWGDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_7_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_7_1/TC_7_1_EKDN220024264-1A_H3MWGDSX5_L3_2.fq.gz &&\ 
echo done
