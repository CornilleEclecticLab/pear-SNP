#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.parse_db_annotation.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/15 17:24:45
# @Description:
#    

import datetime
import sys
import textwrap
import os

from s00_config import *


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

def parse_db(line):
    line = line.strip().replace('>', '')
    
    lines = line.strip().split('|')
    row_id = lines[0].strip()
    
    # Osa from MSU
    if row_id.startswith('LOC') \
        or row_id.startswith('ChrSy') \
        or row_id.startswith('ChrUn'):
        
        gene_id = row_id.split()[0].strip()
        symbol = '-'
        description = lines[1].strip()

    # Ath from TAIR
    elif row_id.startswith('AT'):
        gene_id = row_id
        symbol = lines[1].replace("Symbols: ", "").strip()
        if symbol == '':
            symbol = '-'
        description = lines[2].strip()
    
    # Pyr from UniProtKB
    elif "OS=" in line:
        gene_id = row_id
        symbol = '-'
        description = ' '.join(lines[1].strip().split('OS=')[0].split(' ')[1:])

    else:
        raise ValueError(f'Unrecognized gene ID: >{line}')    
    
    return (gene_id, symbol, description)

for sp, file in ['Ath', Ath_fa], ['Osa', Osa_fa], ['Pyr', Pyr_fa]:
    print(f'Parsing {sp} annotation {file}...')
    with open(os.path.join(db_dir, file), 'r') as fi, \
         open(os.path.join(s03_dir, sp+'.annotation.tsv'), 'w') as fo:
        parses = [parse_db(l) for l in fi.readlines() if l.startswith('>')]
        for gene_id, symbol, description in parses:
            fo.write(f'{gene_id}\t{symbol}\t{description}\n')


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))