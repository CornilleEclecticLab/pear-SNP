#!/bin/bash
#SBATCH -J ch132.s01.raw_qc
#SBATCH -o ch132.s01.raw_qc.out
#SBATCH -e ch132.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/ch132
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch132 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch132/ch132_EKDN220024212-1A_H3MGWDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/ch132 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch132/ch132_EKDN220024212-1A_H3MGWDSX5_L2_2.fq.gz &&\ 
echo done
