#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493995 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493995 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691620/SRR13493995/SRR13493995_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493995 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691620/SRR13493995/SRR13493995_2.fastq.gz &&\
