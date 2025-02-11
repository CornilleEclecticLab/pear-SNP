#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s05.parse_MCScanX.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/10 16:40:20
# @Description:
#    

# Update: v1.1.0 2025-02-11 10:51:42
# - Fix the bug that handle the intersection and union of collinearity pairs.

import datetime
import sys
import textwrap
import os
import glob
from s00_config import *

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.1.0"
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

collinearity_files = glob.glob(os.path.join(s04_dir, "*.collinearity"))

if len(collinearity_files) != 2:
    raise ValueError(f"Unexpected number of collinearity files: {len(collinearity_files)}")


common_collinearity = {}
fo_intersection = open(os.path.join(sub_output_dir, "intersection.collinearity.txt"), 'w', encoding='utf-8')
fo_union = open(os.path.join(sub_output_dir, "union.collinearity.txt"), 'w', encoding='utf-8')
for file in collinearity_files:
    with open(file, 'r') as f:
        inter_switch = None
        basename = os.path.basename(file)
        prefix = os.path.splitext(basename)[0]
        with open(os.path.join(sub_output_dir, f"{prefix}.parsed.collinearity.txt"), 'w', encoding='utf-8') as fo:
            fo_records = {}
            for line in f:
                if line.startswith("###") or line.startswith("# "): # Skip the split lines, paramaters and statistics.
                    continue
                elif line.startswith("## Alignment"):
                    genomes=line.strip().split()[-2]
                    if "&" not in genomes:
                        raise ValueError(f"Unexpected genome names pair: {genomes}")

                    chr1, chr2 = genomes.split("&")
                    genome1, genome2 = chr1[:2], chr2[:2]

                    if genome1 == genome2:
                        inter_switch = False
                    elif genome1 != genome2:
                        inter_switch = True

                elif inter_switch is True:
                    lines=line.strip().split()
                    gene_pairs = f"{lines[-3]}\t{lines[-2]}"
                    if gene_pairs not in fo_records:
                        fo_records[gene_pairs] = True
                        fo.write(f"{gene_pairs}\n")

                        if gene_pairs not in common_collinearity:
                            common_collinearity[gene_pairs] = True
                            fo_union.write(f"{gene_pairs}\n")
                        elif gene_pairs in common_collinearity:
                            fo_intersection.write(f"{gene_pairs}\n")



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))