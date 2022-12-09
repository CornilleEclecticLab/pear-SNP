#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493742 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493742 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691468/SRR13493742/SRR13493742_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR13493742 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA563813_Zhang_2021/SAMN12691468/SRR13493742/SRR13493742_2.fastq.gz &&\
