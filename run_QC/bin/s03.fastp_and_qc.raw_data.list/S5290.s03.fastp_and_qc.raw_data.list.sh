#!/bin/bash
#SBATCH -J S5290.s03.fastp_and_qc.raw_data.list
#SBATCH -o S5290.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e S5290.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5290/S5290_EKDN220024240-1A_H3MGWDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5290/S5290_EKDN220024240-1A_H3MGWDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/clean.S5290_EKDN220024240-1A_H3MGWDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/clean.S5290_EKDN220024240-1A_H3MGWDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/fastp.S5290.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/fastp.S5290.json -R 'S5290' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/clean.S5290_EKDN220024240-1A_H3MGWDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5290/clean.S5290_EKDN220024240-1A_H3MGWDSX5_L3_2.fq.gz &&\
echo done

