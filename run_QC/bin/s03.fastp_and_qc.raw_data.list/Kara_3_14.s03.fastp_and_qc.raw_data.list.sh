#!/bin/bash
#SBATCH -J Kara_3_14.s03.fastp_and_qc.raw_data.list
#SBATCH -o Kara_3_14.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Kara_3_14.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_3_14/Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_3_14/Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/fastp.Kara_3_14.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/fastp.Kara_3_14.json -R 'Kara_3_14' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H33GNDSX5_L4_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_3_14/Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_3_14/Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/fastp.Kara_3_14.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/fastp.Kara_3_14.json -R 'Kara_3_14' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_3_14/clean.Kara_3_14_EKDN220024220-1A_H3N2VDSX5_L2_2.fq.gz &&\
echo done

