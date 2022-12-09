#!/bin/bash
#SBATCH -J P1908.s03.fastp_and_qc.raw_data.list
#SBATCH -o P1908.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e P1908.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H33GNDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H33GNDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H33GNDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H33GNDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/fastp.P1908.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/fastp.P1908.json -R 'P1908' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H33GNDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H33GNDSX5_L2_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/P1908/P1908_EKDN220024266-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/fastp.P1908.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/fastp.P1908.json -R 'P1908' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/P1908/clean.P1908_EKDN220024266-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

