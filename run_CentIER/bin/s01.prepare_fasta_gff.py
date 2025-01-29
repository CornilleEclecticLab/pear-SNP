#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.prepare_fasta_gff.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/01/23 11:20:24
# @Description:
#
# need Python >= 3.8

import datetime
import sys
import os
import glob

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
os.makedirs(sub_script_dir, exist_ok=True)

# Collect fasta files
with open(os.path.join(sub_script_dir, script_basename+".py"), 'w') as f:
    f.write(f"""#!/usr/bin/env python3
from Bio import SeqIO
""")
    fastas = glob.glob(os.path.join(input_dir, '*.fasta'))
    for fa in fastas:
        basename = os.path.basename(fa)
        prefix = os.path.splitext(basename)[0]
        chr_list_file = os.path.join(input_dir, prefix+'.chr_list.txt')
        chr_map_file = os.path.join(input_dir, prefix+'.chrID_map.txt')
        gff_file = os.path.join(input_dir, prefix+'.chr_list.gff')
        
        output_fasta = os.path.join(sub_output_dir, prefix+'.chr_list.fasta')
        output_gff = os.path.join(sub_output_dir, prefix+'.chr_list.gff')
        
        f.write(f"""
# Parse chrID map file
chr_map = dict()
try:
    with open('{chr_map_file}', 'r') as f:
        for line in f:
            old_id, new_id = line.strip().split()
            chr_map[old_id] = new_id
except FileNotFoundError:
    chr_map = dict()

# Parse chr list file
with open('{chr_list_file}', 'r') as f:
    chr_list = f.read().strip().split('\\n')

# Update fasta file
fasta_file = SeqIO.parse('{fa}', 'fasta')
with open('{output_fasta}', 'w') as f:
    for record in fasta_file:
        if record.id in chr_list:
            if record.id in chr_map:
                record.id = chr_map[record.id]
            record.description = ''
            SeqIO.write(record, f, 'fasta')

# Update gff file
with open('{gff_file}', 'r') as fi, open('{output_gff}', 'w') as fo:
    for line in fi:
        if line.strip()=='':
            continue
        if line.startswith('#'):
            fo.write(line)
        else:
            chr= line.split('\\t')[0]
            if chr in chr_list:
                if chr in chr_map:
                    chr = chr_map[chr]
                fo.write('\\t'.join([chr]+line.split('\\t')[1:]))
            else:
                continue

""")

with open(os.path.join(sub_script_dir, script_basename+'.sh'), 'w') as f:
        f.write(f"""#!/usr/bin/bash
#SBATCH -J {script_basename}.sh
#SBATCH -o {script_basename}.%J.out
#SBATCH -e {script_basename}.%J.err
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 00:03:00

module load python

python3 {os.path.join(sub_script_dir, script_basename+'.py')}

""")

os.system(f'cd {sub_script_dir} && sbatch {script_basename}.sh')
print("This script has been submitted to the cluster.")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
