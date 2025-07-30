#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.a1.merge_sweep_bed.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/08 22:07:35
# @Description:
#
# v2.0.0
#   - Changed output filenames.

import datetime
import sys
import textwrap
import os
import subprocess
import glob
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "2.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')


input_bed_file_dir = os.path.join(input_dir,'selection_bed')

# Function to wrap text to 79 characters
def wrap79(text, width=79):
	return textwrap.fill(text, width=width, subsequent_indent=' '*4)

bed_files = glob.glob(os.path.join(input_bed_file_dir, '*.bed'))

bed_file_dic = {}
for bed in bed_files:
    base_name = os.path.basename(bed)
    pars = base_name.split('.')
    method = pars[0]
    pop = pars[2]
    # Skipe the old ouput file of this scritp
    if method.startswith('selection'):
        continue
    # Only keep files that filename starts with Pixy_TajimaD or OmegaPlus
    # if method.startswith('RAiSD') 
    if method.startswith('Pixy_TajimaD') or method.startswith('OmegaPlus'):
        bed_file_dic.setdefault(pop, []).append(bed)

for pop, bed_list in bed_file_dic.items():
    print(pop,':',bed_list)
    base_name = os.path.basename(bed_list[0])
    dir_path = os.path.dirname(bed_list[0])
    pars = base_name.split('.')
    output_name = 'selection.merge_across_methods.' + pop + '.merged.bed'
    output_bed = os.path.join(dir_path, output_name)

    cmd = (
        f'source s00.load_bedtools.sh && '
        f'cat {" ".join(bed_list)} | bedtools sort | bedtools merge > {output_bed}'
    )
    subprocess.run(cmd, shell=True, check=True,  executable='/bin/bash')

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
