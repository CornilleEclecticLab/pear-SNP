#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : test_count.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/10/24 20:27:30
# @Description:
#    

import datetime
import sys
import textwrap
import os
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


with open("../output/pi.txt", "r") as f:
    dict_diff = {}
    dict_comparison = {}
    count_diff = {}
    count_comparison = {}
    for line in f:
        if line.startswith("pop"):
            continue
        ls = line.split()
        pop = ls[0]
        diff = ls[6]
        comparison = ls[7]
        if diff == "NA" or comparison == "NA":
            continue
        count_comparison[pop] = count_comparison.get(pop, 0) + 1
        count_diff[pop] = count_diff.get(pop, 0) + 1
        diff = int(diff)
        comparison = int(comparison)
        dict_diff[pop] = dict_diff.get(pop, 0) + diff
        dict_comparison[pop] = dict_comparison.get(pop, 0) + comparison
    print(dict_diff)
    print(dict_comparison)

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))