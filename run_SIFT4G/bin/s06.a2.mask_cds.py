#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.a2.mask_cds.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/09 15:08:17
# @Description:
#    

import datetime
import sys
import textwrap
import os
import subprocess
from s06_config import *
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
    
    output_masked_cds = cds_masked_bed_dic[ref]

    cmd = (f'. ./s00.load_bedtools.sh && '
            f"""awk 'BEGIN{{OFS = "\\t"}} $3 == "CDS" {{print $1, $4-1, $5}}' {gff_file} | """ 
            f'bedtools intersect -a - -b {mask_pass_bed} -wa | '
            f'bedtools sort |'
            f'bedtools merge > {output_masked_cds}'
    )
    
    print(f'Processing {ref}...')
    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')




end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))