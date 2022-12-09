#!/bin/bash
#SBATCH -J P_dav_V7_1.s03.fastp_and_qc.raw_data.list
#SBATCH -o P_dav_V7_1.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e P_dav_V7_1.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P_dav_V7_1/P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P_dav_V7_1/P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/clean.P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/clean.P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/fastp.P_dav_V7_1.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/fastp.P_dav_V7_1.json -R 'P_dav_V7_1' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/clean.P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P_dav_V7_1/clean.P_dav_V7_1_EKDN220024224-1A_H3MGWDSX5_L2_2.fq.gz &&\
echo done

