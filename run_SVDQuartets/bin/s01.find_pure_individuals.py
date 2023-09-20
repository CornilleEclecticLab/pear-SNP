#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s01.find_pure_individuals.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/15 11:35:21
# @Description:
#     

import datetime
import sys
import textwrap
import os
from types import new_class
from warnings import warn

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


# Constants
version = "1.0.0"
threshold = 1.0 #0.8   # The threshold value that to distinguish the admixture


# File paths
script_basename = 's01.find_pure_individuals'
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
input_file_Q = os.path.join(input_dir,'test.Q')                # Please change this path to your actual input file
input_file_fam = os.path.join(input_dir,'test.fam') 
input_file_id_map = os.path.join(input_dir,'test.id_map')
output_dir = os.path.join(work_dir,'output',script_basename)
output_file = os.path.join(output_dir,'s01.non_admix')
sub_script_dir = os.path.join(work_dir,'bin',script_basename)

os.system(f'mkdir -p {output_dir}')




def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Create lists and dictionaries
ids = []
u_ids = []
m_ids = []
uni_ids = {}

# Read ids from .fam and id_map, 
with open(input_file_fam, 'r') as fi:
    ids = [line.strip().split()[1] for line in fi.readlines()]

with open(input_file_id_map, 'r') as fi:
    lines = fi.read().strip().split('\n')
    m_ids = [i.split()[0] for i in lines]
    u_ids = [i.split()[1] for i in lines]

# Check ids in fam and id_map
if len(ids) != len(u_ids) or len(ids) != len(m_ids):
    warn(wrap(f'ERROR: Incompatible two input files'))
    exit
for id, m_id in zip(ids,m_ids):
    if id != m_id:
        warn(wrap(f'WARNING: Incompatible two input files'))

# Create a dictionary mapping u_ids to ids
uni_ids = {id:u_id for u_id, id in zip(u_ids,m_ids)}


# Find the non-admixture individuals
fo = open(output_file,'w') 
fo2 = open(output_file+'.id_map.txt','w')
with open(input_file_Q, 'r') as fi:
    n = 0
    for line in fi:
        line = line.strip()
        if line == '':
            continue
        admix = True
        n+=1
        line = line.strip()
        Qs = line.split()
        for value in Qs:
            if float(value) >= threshold:
                admix = False
                id = ids[n]
                fo.write(f"{id}\n")
                fo2.write(f"{id}\t{uni_ids[id]}\n")
                break
fo.close()


if n != len(ids):
    print(n, len(ids)-1)
    warn(wrap('WARNING: Incompatible two input files'))
    exit


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))