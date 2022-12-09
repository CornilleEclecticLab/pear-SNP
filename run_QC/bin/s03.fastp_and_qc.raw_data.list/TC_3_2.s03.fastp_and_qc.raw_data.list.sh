#!/bin/bash
#SBATCH -J TC_3_2.s03.fastp_and_qc.raw_data.list
#SBATCH -o TC_3_2.s03.fastp_and_qc.raw_data.list.out
#SBATCH -e TC_3_2.s03.fastp_and_qc.raw_data.list.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2
fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_3_2/TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_3_2/TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/fastp.TC_3_2.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/fastp.TC_3_2.json -R 'TC_3_2' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H35GFDSX5_L1_2.fq.gz &&\
echo done

fastp -f 2 -l 50 \
-i /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_3_2/TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_1.fq.gz -I /work/zruilin/Peach/QC/input/peach_Illumina_novogene_25082022/merged_row_data/TC_3_2/TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_2.fq.gz \
-o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_1.fq.gz -O /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_2.fq.gz \
-h /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/fastp.TC_3_2.html -j /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/fastp.TC_3_2.json -R 'TC_3_2' &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_1.fq.gz &&\
fastqc -t 2 -f fastq -o /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2 /work/zruilin/Peach/QC/output/s03.fastp_and_qc.raw_data.list/TC_3_2/clean.TC_3_2_EKDN220024262-1A_H3MWGDSX5_L3_2.fq.gz &&\
echo done

