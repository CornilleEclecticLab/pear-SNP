#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318832 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318832 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863057/SRR14318832/SRR14318832_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318832 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863057/SRR14318832/SRR14318832_2.fastq.gz &&\
