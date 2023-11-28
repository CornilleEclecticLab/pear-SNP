#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s02.format_pca_data_for_plot.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/22 16:27:47
# @Description:
#     

import datetime
import sys
import textwrap
import os

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


# Global variables
version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(work_dir,'bin',script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


if len(sys.argv) != 4:
    print(wrap(f'''
        Usage: {sys.argv[0]} <id_map> <color_panel> <eigenvec>
        '''))
    print(wrap(f'''
        Example: {sys.argv[0]} example.id_map.txt example.color_panel.txt exmaple.eigenvec
        '''))
    sys.exit(1)


with open(sys.argv[1],'r') as fi:
    id_map = dict()
    lines = fi.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            line = line.split()
            id_map[line[0]] = line[1]


with open(sys.argv[2],'r') as fi:
    color_panel = dict()
    lines = fi.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            line = line.split()
            color_panel[line[0]] = line[1]


with open(sys.argv[3],'r') as fi:
    pca_data = dict()
    lines = fi.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            line = line.replace(' ','\t')
            line = line.split()
            id = line[1]
            pca_data[id] = line
pcs = '\tPC'.join([str(n) for n in range(1,len(line)-1)])
head =f"ID\tGroup\tColor\tPC{pcs}\n"


if len(id_map) != len(pca_data) or len(pca_data) != len(color_panel):
    print('Error: The number of samples in id_map, color_panel and pca_data is not equal!')
    sys.exit(1)


with open(os.path.basename(sys.argv[3])+'.input_pca_lot.data','w') as fo:
    fo.write(head)
    for id, dt in pca_data.items():
        dts = '\t'.join(dt[2:])
        fo.write(f"{id}\t{id_map[id]}\t\"{color_panel[id]}\"\t{dts}\n")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))