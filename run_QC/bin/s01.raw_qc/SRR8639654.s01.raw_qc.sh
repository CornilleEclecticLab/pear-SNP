#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639654 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639654 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA524076/SAMN10995774/SRR8639654/SRR8639654_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR8639654 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA524076/SAMN10995774/SRR8639654/SRR8639654_2.fastq.gz &&\
