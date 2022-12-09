#!/bin/bash
#SBATCH -J S7190.s03.fastp_and_qc.raw_data.list
#SBATCH -o S7190.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e S7190.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/fastp.S7190.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/fastp.S7190.json -R 'S7190' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S7190/S7190_EKDN220024256-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/fastp.S7190.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/fastp.S7190.json -R 'S7190' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S7190/clean.S7190_EKDN220024256-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

