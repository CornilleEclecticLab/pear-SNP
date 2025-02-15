#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s05.mask_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025-01-23 16:18:51
# @Description:
#    Mask the outliers from RAiSD output.


import datetime
import sys
import textwrap
import os
import glob
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')

# The cutoff for overlap with the mask pass region
percent_cutoff = 0.8

# The method
method = 'RAiSD'

version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')
sub_output_dir = os.path.join(work_dir, 'output', script_basename)
sub_script_dir = os.path.join(bin_dir, script_basename)
s04_dir = os.path.join(output_dir, 's04.catch_outliers')

# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)
os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Load populations
populations = {}
with open(os.path.join(input_dir, individual_pop_map_filename)) as f:
    for line in f:
        ind, pop = line.strip().split()
        populations[pop] = populations.get(pop, []) + [ind]


# Load chromosomes
chromosomes = []
with open(os.path.join(input_dir, chromosomes_list_filename)) as f:
    for line in f:
        chromosomes.append(line.strip())
    chromosomes = sorted(chromosomes)   # Sort chromosomes


# Load cutoffs
cutoff_pop = {}
for pop in populations:
    cutoff_file = os.path.join(
        input_dir, cutoff_files_dir, f'{pop}.{cutoff_files_suffix}')
    with open(cutoff_file) as f:
        lines = f.readlines()
        header = lines[0].strip()
        cutoff = lines[-1].strip()
        if header.startswith('#') or header.startswith(pop):
            cutoff_pop[pop] = float(cutoff)
        else:
            raise ValueError(f'Unexpected format in {cutoff_file}')


# Set the gff file path
gff_file = os.path.join(input_dir, gff_filename)

# Collect the outliers files
outliers_files = glob.glob(os.path.join(s04_dir, '*outliers.*.bed'))

# Deal with each outliers file
for file in outliers_files:
    basename = os.path.basename(file)
    if basename.endswith('merged.bed'):
        continue    
    
    print(wrap79(f'Processing {file}'))

    prefix = os.path.splitext(basename)[0]
    pop = basename.split('.')[2]
    
    output_bed = os.path.join(sub_output_dir, f'{prefix}.filtered.sorted.merged.bed')
    output_prefix = os.path.join(sub_output_dir, f'{prefix}')
    
    # Write a sub script
    sub_prefix = f's05.{prefix}.filter'
    sub_script = f'{sub_prefix}.sh'
    with open(os.path.join(sub_script_dir, sub_script ), 'w') as f:
        f.write(f"""#!/bin/bash
#SBATCH -J {sub_prefix}
#SBATCH -o {sub_prefix}.%j.out
#SBATCH -e {sub_prefix}.%j.err

{load_bedtools}
percent_cutoff={percent_cutoff}

bedtools intersect \\
    -a {file} \\
    -b {os.path.join(input_dir, mask_pass_filename)} \\
    -wao \\
| awk -v percent_cutoff="$percent_cutoff" '{{if($7 / ($3 - $2) >= percent_cutoff) print $1"\\t"$2"\\t"$3}}' \\
| bedtools sort \\
| bedtools merge \\
> {output_bed}

cat {output_bed} \\
| awk '{{len=$3-$2;sum+=len}}END{{print sum}}' \\


# Get candidate genes
bedtools intersect \\
    -a {gff_file} \\
    -b {output_bed} \\
    -wa \\
| awk '{{$2="{method}"; print}}' \\
> {output_prefix}.genes.gff

cat {output_prefix}.genes.gff \\
| awk -F '[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\\t{method}"}}' \\
> {output_prefix}.genes.txt

if [[ -s "{output_prefix}.genes.txt" ]]; then
    sort -u {output_prefix}.genes.txt > {output_prefix}.genes.txt.tmp && \\
    mv {output_prefix}.genes.txt.tmp {output_prefix}.genes.txt
fi
""")
        
    os.system(f"cd {sub_script_dir} && sbatch {sub_script}")
    print(wrap79(f"The following script has been submitted to slurm:"))
    print(wrap79(sub_script))
    print('')


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
