#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.x1.mask_gene.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/09 15:08:17
# @Description:
#    Mask specific (e.g. low-mappability) regions in CDS from GFF files using BEDTools.

#  origianl file name ../../run_SIFT4G/bin/s06.a2.mask_cds.py

import datetime
import sys
import textwrap
import os
import subprocess
from s01_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')


for ref in ref_genomes:
    gff_file = gff_file_dic[ref]
    mask_pass_bed = mask_pass_bed_dic[ref]
    
    gene_bed = gene_bed_dic[ref]
    output_masked_gene = gene_masked_bed_dic[ref]

    cmd = (f"module purge && module load bedtools && "
           f"""awk 'BEGIN{{OFS = "\\t"}} $3 == "gene" {{print $1, $4-1, $5}}' {gff_file} > {gene_bed} && """ 
            f'bedtools intersect -a {gene_bed} -b {mask_pass_bed} -wa | '
            f'bedtools sort |'
            f'bedtools merge > {output_masked_gene}'
    )
    
    print(f'Processing {ref} to get bed file of masked gene...')
    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
    
    
    gene_gff = gene_gff_dic[ref]
    output_masked_gene_gff = gene_masked_gff_dic[ref]
    cmd2 = (f"module purge && module load bedtools && "
           f"""awk 'BEGIN{{OFS = "\\t"}} $3 == "gene" {{print}}' {gff_file} > {gene_gff} && """ 
            f'bedtools intersect -a {gene_gff} -b {mask_pass_bed} -wa | '
            f'bedtools sort > {output_masked_gene_gff}'
    )
    
    print(f'Processing {ref} to get gff file of masked gene...')
    subprocess.run(cmd2, shell=True, check=True, executable='/bin/bash')
    





end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))