#!/bin/bash
#SBATCH -J TC_8_1.s03.fastp_and_qc.raw_data.list
#SBATCH -o TC_8_1.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e TC_8_1.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_8_1/TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_8_1/TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/clean.TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/clean.TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/fastp.TC_8_1.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/fastp.TC_8_1.json -R 'TC_8_1' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/clean.TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_8_1/clean.TC_8_1_EKDN220024265-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

