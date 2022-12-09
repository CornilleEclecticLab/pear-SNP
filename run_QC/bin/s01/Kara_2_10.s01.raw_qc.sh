#!/bin/bash
#SBATCH -J Kara_2_10.s01.raw_qc
#SBATCH -o Kara_2_10.s01.raw_qc.out
#SBATCH -e Kara_2_10.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_2_10
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_2_10 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_2_10/Kara_2_10_EKDN220024219-1A_H3N2VDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_2_10 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_2_10/Kara_2_10_EKDN220024219-1A_H3N2VDSX5_L2_2.fq.gz &&\ 
echo done
