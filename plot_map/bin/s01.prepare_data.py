#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s01.prepare_data.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/22 11:02:02
# @Description:
#     

import datetime
import sys
import textwrap
import os
import numpy as np
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
input_chose_id_map = os.path.join(input_dir,'test.txt')
input_locate = os.path.join(input_dir,'test.locate.txt')
input_Ne = os.path.join(input_dir,'test.ne.txt')
output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(work_dir,'bin',script_basename)


# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)


# Function to wrap text to 79 characters
def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Read the id and locate
with open(input_locate,'r') as f:
    dic_locate = {}
    lines = f.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            id = line.split()[0]
            locate = line.split()[1:]
            dic_locate[id] = locate


# Read the id and Q value 
with open(input_chose_id_map,'r') as f:
    chose_IDs = []
    dic_sp_Q = {}
    lines = f.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        parts = line.split()
        id = parts[0]
        sp = parts[1]
        uni_id = parts[2]
        Qs = np.array([float(q) for q in parts[3:]])

        dic_sp_Q[sp].setdefault(sp, value=[]).append(Qs)


# Calculate the average Qs for each species
dic_avg_Q = {sp: np.mean(np.array(Qs), axis=0).tolist() for sp, Qs in dic_sp_Q.items()}


# Check the sum of Qs if equal to 1.0
for Qs in dic_avg_Q.values():
    if np.sum(Qs) != 1.0:
        print('ERROR: The sum of Qs is not 1.0, but {}'.format(np.sum(Qs)))
        sys.exit(1)



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))