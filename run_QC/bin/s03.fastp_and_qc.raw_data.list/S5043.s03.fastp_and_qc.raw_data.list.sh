#!/bin/bash
#SBATCH -J S5043.s03.fastp_and_qc.raw_data.list
#SBATCH -o S5043.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e S5043.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5043/S5043_EKDN220024238-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5043/S5043_EKDN220024238-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/fastp.S5043.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/fastp.S5043.json -R 'S5043' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5043/S5043_EKDN220024238-1A_H3NGMDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/S5043/S5043_EKDN220024238-1A_H3NGMDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H3NGMDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H3NGMDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/fastp.S5043.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/fastp.S5043.json -R 'S5043' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H3NGMDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/S5043/clean.S5043_EKDN220024238-1A_H3NGMDSX5_L3_2.fq.gz &&\
echo done

