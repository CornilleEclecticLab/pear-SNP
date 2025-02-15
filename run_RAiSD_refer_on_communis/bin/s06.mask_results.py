#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.mask_results.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025-01-23 16:34:10
# @Description:
#    Mask the results from RAiSD output.


import datetime
import sys
import textwrap
import os
import glob
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')
sub_output_dir = os.path.join(work_dir, 'output', script_basename)
sub_script_dir = os.path.join(bin_dir, script_basename)


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


# The cutoff for overlap with the mask pass region
percent_cutoff = 0.8

# The method
method = 'RAiSD'

# Set the gff file path
gff_file = os.path.join(input_dir, gff_filename)

# Collect the results files
s03_dir = os.path.join(output_dir, 's03.generate_run_'+method)
result_files = glob.glob(os.path.join(s03_dir, method+'_Report.*'))

# Deal with each results file
half_grid_window_size = grid_window_size // 2
for file in result_files:
    basename = os.path.basename(file)

    print(wrap79(f'Processing {file}'))

    prefix = basename
    pop = basename.split('.')[1]
    chr = basename.split('.')[2]
    output_result = os.path.join(sub_output_dir, f'{basename}.filtered')
    temp_result_bed = os.path.join(sub_output_dir, f'temp.{basename}.bed')
    fo_temp = open(temp_result_bed, 'w')
    with open(file,'r') as fi:
        for line in fi:
            line = line.strip()
            if line.startswith('//') or line == '':
                continue
            pos, start, end, var, sfs, ld, mu = line.split()
            
            pos = int(pos)
            if start == '0' and end == '0':
                start = max(pos - half_grid_window_size, 0)
                end = pos + half_grid_window_size
            else:
                start = int(start)-1
                end = int(end)
            
            fo_temp.write(
                f'{chr}\t{start}\t{end}\t{pos}\t{var}\t{sfs}\t{ld}\t{mu}\n')
    fo_temp.close()
    # Write a sub script
    sub_prefix = f's06.{prefix}.filter'
    sub_script = f'{sub_prefix}.filter.sh'
    with open(os.path.join(sub_script_dir, sub_script ), 'w') as f:
        f.write(f"""#!/bin/bash
#SBATCH -J {sub_prefix}
#SBATCH -o {sub_prefix}.%j.out
#SBATCH -e {sub_prefix}.%j.err

{load_bedtools}
percent_cutoff={percent_cutoff}

bedtools intersect \\
    -a {temp_result_bed} \\
    -b {os.path.join(input_dir, mask_pass_filename)} \\
    -wao \\
| awk -v percent_cutoff="$percent_cutoff" '{{if($12 / ($3 - $2) >= percent_cutoff) print $4"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$7"\\t"$8}}' \\
> {output_result}

cat {output_result} \\
| awk '{{len=$3-$2;sum+=len}}END{{print sum}}' \\
""")
        
    os.system(f"cd {sub_script_dir} && sbatch {sub_script}")
    print(wrap79(f"The following script has been submitted to slurm:"))
    print(wrap79(sub_script))
    print('')


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
