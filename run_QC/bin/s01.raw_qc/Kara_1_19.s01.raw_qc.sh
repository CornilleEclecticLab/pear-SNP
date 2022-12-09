#!/bin/bash
#SBATCH -J Kara_1_19.s01.raw_qc
#SBATCH -o Kara_1_19.s01.raw_qc.out
#SBATCH -e Kara_1_19.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_1_19
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_1_19 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_1_19 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_1_19 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Kara_1_19 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_2.fq.gz &&\ 
echo done
