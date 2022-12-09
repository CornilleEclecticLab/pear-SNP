#!/bin/bash
#SBATCH -J TC_8_1.s01.raw_qc
#SBATCH -o TC_8_1.s01.raw_qc.out
#SBATCH -e TC_8_1.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/TC_8_1
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_8_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_8_1/TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_8_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_8_1/TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_2.fq.gz &&\ 
echo done
