#!/bin/bash
#SBATCH -J Clone_5.s03.fastp_and_qc.raw_data.list
#SBATCH -o Clone_5.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e Clone_5.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_5/Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/Clone_5/Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/clean.Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/clean.Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/fastp.Clone_5.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/fastp.Clone_5.json -R 'Clone_5' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/clean.Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/Clone_5/clean.Clone_5_EKDN220024217-1A_H3N2VDSX5_L2_2.fq.gz &&\
echo done

