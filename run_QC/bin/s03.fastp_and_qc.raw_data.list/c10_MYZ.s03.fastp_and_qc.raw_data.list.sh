#!/bin/bash
#SBATCH -J c10_MYZ.s03.fastp_and_qc.raw_data.list
#SBATCH -o c10_MYZ.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e c10_MYZ.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir -p /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/c10_MYZ/c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/c10_MYZ/c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/fastp.c10_MYZ.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/fastp.c10_MYZ.json -R 'c10_MYZ' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/c10_MYZ/c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/c10_MYZ/c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/fastp.c10_MYZ.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/fastp.c10_MYZ.json -R 'c10_MYZ' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/c10_MYZ/clean.c10_MYZ_EKDN220024267-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

