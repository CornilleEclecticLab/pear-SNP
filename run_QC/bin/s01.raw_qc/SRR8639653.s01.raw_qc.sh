#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639653 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639653 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA524076/SAMN10995775/SRR8639653/SRR8639653_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639653 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA524076/SAMN10995775/SRR8639653/SRR8639653_2.fastq.gz &&\
