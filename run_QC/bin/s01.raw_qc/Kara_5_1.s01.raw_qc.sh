#!/bin/bash
#SBATCH -J Kara_5_1.s01.raw_qc
#SBATCH -o Kara_5_1.s01.raw_qc.out
#SBATCH -e Kara_5_1.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_5_1
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_5_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_5_1/Kara_5_1_EKDN220024223-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_5_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_5_1/Kara_5_1_EKDN220024223-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_5_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_5_1/Kara_5_1_EKDN220024223-1A_H3N2VDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_5_1 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_5_1/Kara_5_1_EKDN220024223-1A_H3N2VDSX5_L2_2.fq.gz &&\ 
echo done
