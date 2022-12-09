#!/bin/bash
#SBATCH -J S3359.s01.raw_qc
#SBATCH -o S3359.s01.raw_qc.out
#SBATCH -e S3359.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/S3359
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S3359 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S3359/S3359_EKDN220024234-1A_H3MGWDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S3359 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S3359/S3359_EKDN220024234-1A_H3MGWDSX5_L3_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S3359 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S3359/S3359_EKDN220024234-1A_H3NGMDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/S3359 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S3359/S3359_EKDN220024234-1A_H3NGMDSX5_L3_2.fq.gz &&\ 
echo done
