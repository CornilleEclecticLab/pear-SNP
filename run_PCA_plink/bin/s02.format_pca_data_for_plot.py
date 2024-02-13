#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.format_pca_data_for_plot.py
# @Version  : 2.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/22 16:27:47
# @Description:
#     Update version to 2.0.0 2024-02-07 11:33:02
#       1. Using argparse to parse command line arguments
#       2. Allow id_map contains full dataset
#       3. Output file contains more information

import datetime
import sys
import textwrap
import os
import argparse

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


# Global variables
version = "2.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output', script_basename)
sub_script_dir = os.path.join(work_dir, 'bin', script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Read arguments
parser = argparse.ArgumentParser(
    description='Format plink PCA result for plot.')

parser.add_argument(
    '-e', '--eigenvec',
    required=True,
    help='.eigenvec file from plink PCA result.')

parser.add_argument(
    '-i', '--id_map',
    required=True,
    help='ID map file. Content: Sample_ID in vcf file <tab> Unified_ID '
         '[<tab> Group]. The Unified ID expected contains the group name, '
         'e.g. "Species_Country_Group122", or give the group name in the '
         'third column.')

parser.add_argument(
    '-c', '--color_panel',
    required=True,
    help='Color panel file. This is a no header file, contains color in '
         'HEX format with #, each line for one sample, the samples order '
         'should follow the plink PCA results.')

args = parser.parse_args()


# Read input files
with open(args.id_map, 'r') as fi:
    id_map = dict()
    id_map2 = dict()
    group_map = dict()
    lines = fi.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            if len(line.split()) == 2:
                a, b = line.split()
                id_map[a] = b
                id_map2[b] = a

                c = b.split('_')[-1][:2]
                group_map[b] = c
                group_map[a] = c

            elif len(line.split()) == 3:
                a, b, c = line.split()
                id_map[a] = b
                id_map2[b] = a
                group_map[a] = c
                group_map[b] = c

with open(args.color_panel, 'r') as fi:
    color_panel = fi.read().strip().splitlines()

with open(args.eigenvec, 'r') as fi:
    pca_data = dict()
    lines = fi.read().strip().splitlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            line = line.replace(' ', '\t')
            line = line.split()
            id = line[1]
            pca_data[id] = line[2:]
pcs = '\tPC'.join([str(n) for n in range(1, len(line)-1)])
##### To update
header = f"ID_vcf\tID\tGroup\tColor\tPC{pcs}\n"


# Check data
for id in pca_data.keys():
    if id not in id_map.keys():
        for id2 in pca_data.keys():
            if id2 not in id_map2.keys():
                raise ValueError(
                    f'Error: {id} and {id2} not match id_map!')
        else:
            id_map = id_map2

        break

if len(pca_data) != len(color_panel):
    print(
        wrap('Error: The number of samples not equal in color_panel and'
             'pca_data!')
          )
    sys.exit(1)


# Output data
with open(args.eigenvec+'.input_pca_lot.data', 'w') as fo:
    fo.write(header)
    for i, (id, pca_dts) in enumerate(pca_data.items()):
        pca_dt = '\t'.join(pca_dts)
        group = group_map[id]
        color = color_panel[i]
        if color == '#888888':
            group = 'Admixed'
        fo.write(f'{id}\t{id_map[id]}\t{group}\t"{color}"\t{pca_dt}\n')


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
