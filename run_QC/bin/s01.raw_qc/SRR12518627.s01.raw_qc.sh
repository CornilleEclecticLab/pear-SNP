#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR12518627 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR12518627 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA607895/SAMN14485717/SRR12518627/SRR12518627_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR12518627 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA607895/SAMN14485717/SRR12518627/SRR12518627_2.fastq.gz &&\
