#!/bin/bash
#SBATCH -J S8297.s03.fastp_and_qc.raw_data.list
#SBATCH -o S8297.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e S8297.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8297/S8297_EKDN220024257-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8297/S8297_EKDN220024257-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/fastp.S8297.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/fastp.S8297.json -R 'S8297' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8297/S8297_EKDN220024257-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S8297/S8297_EKDN220024257-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/fastp.S8297.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/fastp.S8297.json -R 'S8297' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S8297/clean.S8297_EKDN220024257-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

