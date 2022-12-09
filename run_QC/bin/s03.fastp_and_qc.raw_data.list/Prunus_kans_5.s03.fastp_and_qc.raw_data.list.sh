#!/bin/bash
#SBATCH -J Prunus_kans_5.s03.fastp_and_qc.raw_data.list
#SBATCH -o Prunus_kans_5.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Prunus_kans_5.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/fastp.Prunus_kans_5.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/fastp.Prunus_kans_5.json -R 'Prunus_kans_5' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3MGWDSX5_L2_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Prunus_kans_5/Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/fastp.Prunus_kans_5.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/fastp.Prunus_kans_5.json -R 'Prunus_kans_5' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Prunus_kans_5/clean.Prunus_kans_5_EKDN220024229-1A_H3NGMDSX5_L3_2.fq.gz &&\
echo done

