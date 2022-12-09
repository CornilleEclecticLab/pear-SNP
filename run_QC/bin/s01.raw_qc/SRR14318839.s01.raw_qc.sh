#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318839 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318839 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863042/SRR14318839/SRR14318839_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318839 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863042/SRR14318839/SRR14318839_2.fastq.gz &&\
