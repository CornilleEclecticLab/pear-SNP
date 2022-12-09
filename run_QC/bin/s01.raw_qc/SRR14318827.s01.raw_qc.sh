#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318827 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318827 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863044/SRR14318827/SRR14318827_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR14318827 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA724820/SAMN18863044/SRR14318827/SRR14318827_2.fastq.gz &&\
