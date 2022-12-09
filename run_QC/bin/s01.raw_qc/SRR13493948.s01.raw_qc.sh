#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493948 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493948 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691662/SRR13493948/SRR13493948_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493948 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691662/SRR13493948/SRR13493948_2.fastq.gz &&\
