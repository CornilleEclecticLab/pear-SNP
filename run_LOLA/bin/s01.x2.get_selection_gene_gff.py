#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.x2.get_sweep_gff.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/08/11 13:59:55
# @Description:
#    

import datetime
import sys
import textwrap
import os
import subprocess
import glob
from s01_config import selective_gene_gff_list, summary_lists, wild_pops
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')


# Function to wrap text to 79 characters
def wrap79(text, width=79):
	return textwrap.fill(text, width=width, subsequent_indent=' '*4)


for ref, ls in summary_lists.items():
    with open(os.path.join(input_dir, ls), 'r') as f:
        for line in f:
            path = line.strip()
            path = path.replace('.txt', '.gff')
            path = os.path.realpath(path)
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found: {path}")
            filename = os.path.basename(path)
            pop = filename.split('.')[2]
            print(pop)

            if pop in wild_pops:
                # print(f'Pop {pop} is a wild population, skipping...')
                continue # skip wild populations
            else:
                print(f'Pop {pop} -> {ref}')
                selective_gene_gff_list.setdefault(ref, []).append(path)


for ref, gffs in selective_gene_gff_list.items():
    print(f'Processing {ref}...')
    out_gff_lines = []
    for gff in gffs:
        print(f'Processing GFF file: {gff}')
        in_gff = open(gff, 'r')
        for line in in_gff:
            if ' ' in line:
                raise ValueError("Invalid GFF format, as space exist")
            out_gff_lines.append(line)

    out_gff_lines = set(out_gff_lines)
    output_gff_file = os.path.join(input_dir, 'gff', f'{ref}.selection.gff')
    with open(output_gff_file, 'w') as out_file:
        out_file.write(''.join(out_gff_lines))
    output_bed_file = os.path.join(input_dir, 'bed', f'{ref}.selection.bed')
    with open(output_bed_file, 'w') as out_file:
        for line in out_gff_lines:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                raise ValueError("Invalid GFF format, as fields are missing")
            chrom, start, end, name = fields[0], int(fields[3])-1, int(fields[4]), fields[8]
            out_file.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    cmd = f"module purge && module load bedtools && bedtools sort -i {output_gff_file} > {output_gff_file}.sorted && mv {output_gff_file}.sorted {output_gff_file} && bedtools sort -i {output_bed_file} > {output_bed_file}.sorted && mv {output_bed_file}.sorted {output_bed_file}"
    subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
