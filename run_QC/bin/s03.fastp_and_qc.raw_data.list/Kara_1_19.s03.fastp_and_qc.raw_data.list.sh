#!/bin/bash
#SBATCH -J Kara_1_19.s03.fastp_and_qc.raw_data.list
#SBATCH -o Kara_1_19.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Kara_1_19.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/fastp.Kara_1_19.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/fastp.Kara_1_19.json -R 'Kara_1_19' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3MGWDSX5_L2_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_1_19/Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/fastp.Kara_1_19.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/fastp.Kara_1_19.json -R 'Kara_1_19' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_1_19/clean.Kara_1_19_EKDN220024218-1A_H3NGMDSX5_L4_2.fq.gz &&\
echo done

