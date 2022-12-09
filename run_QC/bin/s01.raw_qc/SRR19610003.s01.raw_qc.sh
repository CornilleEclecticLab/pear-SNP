#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR19610003 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR19610003 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA846875_Qin_2022/SAMN28906012/SRR19610003/SRR19610003_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR19610003 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA846875_Qin_2022/SAMN28906012/SRR19610003/SRR19610003_2.fastq.gz &&\
