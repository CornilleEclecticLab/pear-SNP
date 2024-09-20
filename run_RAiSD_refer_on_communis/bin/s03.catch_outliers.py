#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.catch_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024-09-18 16:43:01
# @Description:
#    Catch outliers in RAiSD output.

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
sub_output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)
os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Load populations
populations = {}
with open(os.path.join(input_dir,individual_pop_map_filename)) as f:
    for line in f:
        ind, pop = line.strip().split()
        populations[pop] = populations.get(pop, []) + [ind]


# Load chromosomes
chromosomes = []
with open(os.path.join(input_dir,chromosomes_list_filename)) as f:
    for line in f:
        chromosomes.append(line.strip())
    chromosomes = sorted(chromosomes)


# Load cutoffs
cutoff_pop = {}
for pop in populations:
    cutoff_file = os.path.join(input_dir,cutoff_files_dir,f'{pop}.{cutoff_files_suffix}')
    with open(cutoff_file) as f:
        lines = f.readlines()
        header = lines[0].strip()
        cutoff = lines[-1].strip()
        if header.startswith('omegas') or header.startswith('#') or header.startswith(pop):
            cutoff_pop[pop] = float(cutoff)
        else:
            raise ValueError(f'Unexpected format in {cutoff_file}')


# Catch outliers from RAiSD by grid output
half_grid_window_size = grid_window_size // 2
pop_outliers_count = {}
for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'RAiSD.outliers.{pop}.bed'), 'w')
    outliers_count_pop = 0
    for chr in chromosomes:
        outliers = []
        # Read the RAiSD output
        raisd_output_file = os.path.join(output_dir,raisd_output_dir,raisd_files_name.format(pop=pop,chr=chr))
        with open(raisd_output_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line == '':
                    continue
                if line.startswith('//') or line.startswith('Position'):
                    continue
                else:
                    try:
                        position, start, end, var, sfs, ld, mu  = line.split()
                    except ValueError:
                        try:
                            position, mu = line.split()
                            start, end = 0, 0
                        except ValueError:
                            raise ValueError(f'Unexpected format in {raisd_output_file}\n{line}')
                    mu = float(mu)
                    position = int(float(position))
                    start = int(float(start))
                    end = int(float(end))
                    if mu >= cutoff_pop[pop]:
                        if start == 0 and end == 0:
                            start = max(position - half_grid_window_size - 1, 0)
                            end = position + half_grid_window_size
                        if (chr, start, end) not in outliers:
                            outliers.append((chr, start, end))
                        else:
                            print(f'Outlier {chr}:{position} {start} {end} already found in {pop}')
        # sort the outliers by start position
        if len(outliers) > 0:
            outliers.sort(key=lambda x: x[1])    # Note: should sort chromosome first, sort -k1,1 -k2,2n
            outliers_count_pop += len(outliers)
            for chr, start, end in outliers:
                output_file.write(f'{chr}\t{start}\t{end}\n')
    if outliers_count_pop == 0:
        print(f'No outlier was found in {pop}')
    output_file.close()

# Merge the outliers
merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge.sh')
with open(merge_script,'w') as f:
    f.write(f"#!/bin/bash\n")
    f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
    f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.out\n\n")
    f.write(f"{load_bedtools}\n\n\n")
    for pop in populations:
        input_file = os.path.join(sub_output_dir, f'RAiSD.outliers.{pop}.bed')
        output_file = os.path.join(sub_output_dir, f'RAiSD.outliers.{pop}.merged.bed')
        f.write(f"bedtools merge -d {grid_window_size} -i {input_file} > {output_file}\n\n")

os.system(f'cd {sub_script_dir} && sbatch {merge_script}')
print(wrap79(f"The following script has been submitted to slurm:"))
print(f'{merge_script}')


# Catch outliers from RAiSD by window output
threshold = 0.001  # Top 1% of the distribution

for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'RAiSD.outliers.by_win.{pop}.bed'), 'w')
    outliers_pop = []
    mu_count = 0
    for chr in chromosomes:
        results = []
        # Read the RAiSD output
        raisd_output_file = os.path.join(
            output_dir, raisd_output_dir, raisd_files_by_win_name.format(pop=pop, chr=chr))
        with open(raisd_output_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line == '':
                    continue
                if line.startswith('//') or line.startswith('Position'):
                    continue
                try:
                    position, start, end, var, sfs, ld, mu = line.split()
                except ValueError:
                    try:
                        position, mu = line.split()
                    except ValueError:
                        raise ValueError(
                            f'Unexpected format in {raisd_output_file}\n{line}')
                mu = float(mu)
                mu_count += 1
                start = int(float(start))
                end = int(float(end))
                if mu > 0:
                    results.append((chr,start, end, mu))
                else:
                    raise ValueError(f'Unexpected mu value in {raisd_output_file}\n{line}')
        # sort the outliers by start position
        if len(results) > 0:
            # Note: should sort chromosome first, sort -k1,1 -k2,2n
            outliers_pop.extend(results)
    outliers_pop.sort(key=lambda x: x[3], reverse=True)
    outliers_pop = outliers_pop[:int(len(outliers_pop)*threshold)]
    outliers_pop.sort(key=lambda x: (x[0], x[1]))

    for chr, start, end, mu in outliers_pop:
        output_file.write(f'{chr}\t{start}\t{end}\n')
    output_file.close()

# Merge the outliers
merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge_by_win.sh')
with open(merge_script, 'w') as f:
    f.write(f"#!/bin/bash\n")
    f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
    f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.out\n\n")
    f.write(f"{load_bedtools}\n\n\n")
    for pop in populations:
        input_file = os.path.join(
            sub_output_dir, f'RAiSD.outliers.by_win.{pop}.bed')
        output_file = os.path.join(
            sub_output_dir, f'RAiSD.outliers.by_win.{pop}.merged.bed')
        f.write(
            f"bedtools merge -d {grid_window_size} -i {input_file} > {output_file}\n\n")

os.system(f'cd {sub_script_dir} && sbatch {merge_script}')
print(wrap79(f"The following script has been submitted to slurm:"))
print(f'{merge_script}')


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
