#!/bin/bash
#SBATCH -J Prunus_kans_5.s01.raw_qc
#SBATCH -o Prunus_kans_5.s01.raw_qc.out
#SBATCH -e Prunus_kans_5.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_5
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_5 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_2.fq.gz &&\ 
echo done
