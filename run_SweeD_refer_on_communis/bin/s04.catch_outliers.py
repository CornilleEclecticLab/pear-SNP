#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.catch_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024-09-18 16:39:51
# @Description:
#    
# Update: v2.0.0 2024-10-01 23:21:07
#    1. Get candidate genes from the outliers.

# Update: v2.1.0 2024-10-08 21:11:00
#    1. Use more cutoffs for outliers.
#    2. Write all valid results to a file.

import datetime
import sys
import textwrap
import os
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "2.1.0"
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


for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'SweeD.outliers.{pop}.bed'), 'w')
    results_pop = []
    for chr in chromosomes:
        outliers = []
        # Read the SweeD output
        sweed_output_file = os.path.join(output_dir,sweed_output_dir,f'SweeD_Report.{pop}.{chr}')
        with open(sweed_output_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line == '':
                    continue
                elif line.startswith('//') or line.startswith('Position'):
                    continue
                else:
                    try:
                        position, likelihood, alpha, start, end = line.split()
                    except ValueError:
                        raise ValueError(f'Unexpected format in {sweed_output_file}\n{line}')
                    
                    clr = float(likelihood)
                    start = int(float(start)) - 1
                    end = int(float(end))
                    
                    if clr > 0:
                        results_pop.append((chr, start, end, clr, position))
                        if clr >= cutoff_pop[pop]:
                            outliers.append((chr, start, end))
            # sort the outliers by start position
            outliers.sort(key=lambda x: x[1])    # Note: should sort chromosome first, sort -k1,1 -k2,2n
            for chr, start, end in outliers:
                output_file.write(f'{chr}\t{start}\t{end}\n')
    output_file.close()
    
    
    # Remove duplicates from results
    results_pop = list(set(results_pop))
    
    # write all valid results to a file
    with open(os.path.join(sub_output_dir, f'SweeD.all.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        results_pop.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{clr}' for (chr, start, end, clr, position) in results_pop]))
    
    
    # Sort the results by clr
    results_pop.sort(key=lambda x: x[3], reverse=True)
    
    # Write the top 30% results
    top_30_results = results_pop[:int(0.3 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'SweeD.top30per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_30_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_30_results]))
    
    # Write the top 20% results
    top_20_results = results_pop[:int(0.2 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'SweeD.top20per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_20_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_20_results]))

    
    top_cut_off = [0.001, 0.01, 200, 500]
    for top in top_cut_off:
        if top < 1:
            tops = max(1,int(float(top * len(results_pop))))
        else:
            tops = int(top)
        tops_outliers = results_pop[:tops]
        tops_outliers.sort(key=lambda x: (x[0], x[1]))
        with open(os.path.join(sub_output_dir, f'SweeD.{top}outliers.{pop}.bed'), 'w') as f:
            for chr, start, end, clr, position in tops_outliers:
                f.write(f'{chr}\t{start}\t{end}\n')
    

# Merge the outliers
def merge_outliers(method, cutoff):
    # cutoff likes '', '0.001', '0.01' or '500'
    global sub_script_dir, script_basename, populations, load_bedtools,input_dir,sub_output_dir,gff_filename
    merge_script = os.path.join(sub_script_dir, f'{script_basename}.{cutoff}merge.sh')
    with open(merge_script,'w') as f:
        f.write(f"#!/bin/bash\n")
        f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
        f.write(f"#SBATCH -o {os.path.basename(merge_script)}.%J.out\n")
        f.write(f"{load_bedtools}\n\n\n")
        for pop in populations:
            input_file = os.path.join(sub_output_dir, f'{method}.{cutoff}outliers.{pop}.bed')
            output_file = os.path.join(sub_output_dir, f'{method}.{cutoff}outliers.{pop}.merged.bed')
            gff_file = os.path.join(input_dir,gff_filename)
            f.write(f"bedtools merge -i {input_file} > {output_file}\n")

            # Get candidate genes
            f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
    | awk '{{$2="{method}"; print}}' \\
    > {output_file}.genes.gff

    cat {output_file}.genes.gff | awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
    > {output_file}.genes.txt
    """)

    os.system(f"cd {sub_script_dir} && sbatch {merge_script}")
    print(wrap79(f"The following script has been submitted to slurm:"))
    print(wrap79(merge_script))

ls = top_cut_off.copy()
ls.append('')
for cutoff in ls:
    merge_outliers('SweeD', cutoff)


# # Merge the outliers
# merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge.sh')
# with open(merge_script,'w') as f:
#     f.write(f"#!/bin/bash\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.out\n")
#     f.write(f"{load_bedtools}\n\n\n")
#     for pop in populations:
#         input_file = os.path.join(sub_output_dir, f'SweeD.outliers.{pop}.bed')
#         output_file = os.path.join(sub_output_dir, f'SweeD.outliers.{pop}.merged.bed')
#         f.write(f"bedtools merge -i {input_file} > {output_file}\n\n")

#         # Get candidate genes
#         gff_file = os.path.join(input_dir, gff_filename)
#         method = "SweeD"
#         f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
# | awk '{{$2="{method}"; print}}' \\
# > {output_file}.genes.gff

# cat {output_file}.genes.gff | awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
# > {output_file}.genes.txt
# """)
# 
# os.system(f'cd {sub_script_dir} && sbatch {merge_script}')
# print(wrap79(f"The following script has been submitted to slurm:"))
# print(wrap79(merge_script))



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))