#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.split_pixy_by_pop.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/07/20 23:29:38
# @Description:
#    

import datetime
import sys
import textwrap
import os
import glob
import threading

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(output_dir,script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)
# os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


tajimaD_dir = os.path.join(
    output_dir, 's01.get_vcf_by_pop_from_allsites_and_cal_TajimaD')
tajimaD_files = glob.glob(os.path.join(tajimaD_dir, '*_tajima_d.txt'))


def process_file(file):
    basename = os.path.basename(file)
    output_name = '.'.join(basename.split('.')[2:])

    with open(file, 'r') as f:
        lines = f.readlines()
        if not lines or not lines[0].startswith('pop'):
            raise ValueError(f"[ERROR] {file} missing header")
            
        header = lines[0].strip()
        headers = header.split()
        pop_index = headers.index('pop')
        tajima_d_index = headers.index('tajima_d')

        pop_lines = {}
        for line in lines[1:]:
            cols = line.strip().split()
            if len(cols) <= max(pop_index, tajima_d_index):
                continue
            if cols[tajima_d_index] == 'NA':
                continue
            pop = cols[pop_index]
            pop_lines.setdefault(pop, []).append(line.strip())

        for pop, content in pop_lines.items():
            out_path = os.path.join(
                sub_output_dir, f"Pixy.{pop}.{output_name}")
            with open(out_path, 'w') as fo:
                fo.write(header + '\n')
                fo.write('\n'.join(content) + '\n')
    
        print(f"Processed {file} into {len(pop_lines)} populations.\n")


threads = []
for file in tajimaD_files:
    t = threading.Thread(target=process_file, args=(file,))
    threads.append(t)
    t.start()

for t in threads:
    t.join()

print("All files processed.")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))