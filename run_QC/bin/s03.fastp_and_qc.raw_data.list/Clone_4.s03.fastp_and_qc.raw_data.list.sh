#!/bin/bash
#SBATCH -J Clone_4.s03.fastp_and_qc.raw_data.list
#SBATCH -o Clone_4.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Clone_4.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/fastp.Clone_4.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/fastp.Clone_4.json -R 'Clone_4' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_4/Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/fastp.Clone_4.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/fastp.Clone_4.json -R 'Clone_4' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_4/clean.Clone_4_EKDN220024216-1A_H3N2VDSX5_L1_2.fq.gz &&\
echo done

