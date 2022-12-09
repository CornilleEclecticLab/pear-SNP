#!/bin/bash
#SBATCH -J Kara_4_4.s03.fastp_and_qc.raw_data.list
#SBATCH -o Kara_4_4.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Kara_4_4.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_4_4/Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_4_4/Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/fastp.Kara_4_4.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/fastp.Kara_4_4.json -R 'Kara_4_4' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_4_4/Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Kara_4_4/Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/fastp.Kara_4_4.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/fastp.Kara_4_4.json -R 'Kara_4_4' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Kara_4_4/clean.Kara_4_4_EKDN220024222-1A_H3N2VDSX5_L2_2.fq.gz &&\
echo done

