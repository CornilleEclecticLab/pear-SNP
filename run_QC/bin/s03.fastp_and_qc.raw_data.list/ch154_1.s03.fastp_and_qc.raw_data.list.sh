#!/bin/bash
#SBATCH -J ch154_1.s03.fastp_and_qc.raw_data.list
#SBATCH -o ch154_1.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e ch154_1.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H33GNDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H33GNDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H33GNDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H33GNDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/fastp.ch154_1.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/fastp.ch154_1.json -R 'ch154_1' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H33GNDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H33GNDSX5_L2_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/ch154_1/ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/fastp.ch154_1.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/fastp.ch154_1.json -R 'ch154_1' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/ch154_1/clean.ch154_1_EKDN220024214-1A_H3N2VDSX5_L1_2.fq.gz &&\
echo done

