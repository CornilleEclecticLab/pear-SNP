#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493979 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493979 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691635/SRR13493979/SRR13493979_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493979 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691635/SRR13493979/SRR13493979_2.fastq.gz &&\
