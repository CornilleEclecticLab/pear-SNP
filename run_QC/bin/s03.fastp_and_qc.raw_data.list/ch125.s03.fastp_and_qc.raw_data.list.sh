#!/bin/bash
#SBATCH -J ch125.s03.fastp_and_qc.raw_data.list
#SBATCH -o ch125.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e ch125.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch125/ch125_EKDN220024211-1A_H3MGWDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch125/ch125_EKDN220024211-1A_H3MGWDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3MGWDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3MGWDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/fastp.ch125.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/fastp.ch125.json -R 'ch125' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3MGWDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3MGWDSX5_L2_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch125/ch125_EKDN220024211-1A_H3TGLDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch125/ch125_EKDN220024211-1A_H3TGLDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3TGLDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3TGLDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/fastp.ch125.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/fastp.ch125.json -R 'ch125' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3TGLDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch125/clean.ch125_EKDN220024211-1A_H3TGLDSX5_L3_2.fq.gz &&\
echo done

