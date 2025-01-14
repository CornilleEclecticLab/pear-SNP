#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.count_hits.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/01/14 17:47:04
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
# os.makedirs(sub_output_dir, exist_ok=True)
# os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Count hits
count_dict = {}
with open(os.path.join(output_dir,"aln.id.txt"),'r') as f:
    for line in f:
        query_id, target_id = line.strip().split()
        count_dict[query_id] = count_dict.get(query_id,{})
        count_dict[query_id][target_id] = count_dict[query_id].get(target_id,0) + 1

# Sort by count
sorted_counts = []
for query_id in count_dict:
    for target_id in count_dict[query_id]:
        sorted_counts.append((query_id,target_id,count_dict[query_id][target_id]))

sorted_counts.sort(key=lambda x: x[2], reverse=True)

# Write to file
with open(os.path.join(output_dir,"aln.id.count.txt"),'w') as f:
    for query_id, target_id, count in sorted_counts:
        f.write(f"{query_id}\t{target_id}\t{count}\n")

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))