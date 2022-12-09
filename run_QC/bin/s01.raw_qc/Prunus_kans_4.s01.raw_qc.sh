#!/bin/bash
#SBATCH -J Prunus_kans_4.s01.raw_qc
#SBATCH -o Prunus_kans_4.s01.raw_qc.out
#SBATCH -e Prunus_kans_4.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_4
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_4/Prunus_kans_4_EKDN220024228-1A_H3N2VDSX5_L2_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/Prunus_kans_4 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_4/Prunus_kans_4_EKDN220024228-1A_H3N2VDSX5_L2_2.fq.gz &&\ 
echo done
