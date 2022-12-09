#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

mkdir -p /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR10556677 
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR10556677 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA592438/SAMN13426700/SRR10556677/SRR10556677_1.fastq.gz &&\
fastqc -t 2 -f -o /data/atipe-workspace/ynie/pear/run_QC/output/s01.raw_qc/SRR10556677 /data/atipe/ynie/catalogue/resequencing_pear/PRJNA592438/SAMN13426700/SRR10556677/SRR10556677_2.fastq.gz &&\
