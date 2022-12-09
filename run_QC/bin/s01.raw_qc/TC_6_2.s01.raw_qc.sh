#!/bin/bash
#SBATCH -J TC_6_2.s01.raw_qc
#SBATCH -o TC_6_2.s01.raw_qc.out
#SBATCH -e TC_6_2.s01.raw_qc.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir /work/zruilin/Peach/QC/output/s01.raw_qc/TC_6_2
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_6_2 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_6_2/TC_6_2_EKDN220024263-1A_H35GFDSX5_L1_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_6_2 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_6_2/TC_6_2_EKDN220024263-1A_H35GFDSX5_L1_2.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_6_2 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_6_2/TC_6_2_EKDN220024263-1A_H3MWGDSX5_L3_1.fq.gz &&\ 
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s01.raw_qc/TC_6_2 /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_6_2/TC_6_2_EKDN220024263-1A_H3MWGDSX5_L3_2.fq.gz &&\ 
echo done
