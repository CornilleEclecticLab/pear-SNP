#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.catch_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/09/12 18:40:50
# @Description:
#    Catch the outliers from OmegaPlus output.

# Update: v1.0.1 2024-09-18 17:20:50
#    1. Resolve the issue when convert a float in string type to int.

# Update: v2.0.0 2024-09-29 21:54:15
#    1. Get candidate genes from the outliers.

# Update: v2.1.0 2024-10-08 15:02:42
#    1. Use more cutoffs for outliers.

# Update: v2.2.0 2024-10-10 10:21:27
#    1. Write all valid results to a file.
#    2. Fix the bug that using string to get boolean value.

import datetime
import sys
import textwrap
import os
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "2.2.0"
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
    chromosomes = sorted(chromosomes)   # Sort chromosomes


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


# Read outputs from OmegaPlus, and catch the outliers based on neutral envlope
for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'OmegaPlus.outliers.{pop}.bed'), 'w')
    results_pop = []
    for chr in chromosomes:
        outliers = []
        # Read the OmegaPlus output
        omega_plus_output_file = os.path.join(output_dir,omega_plus_output_dir,f'OmegaPlus_Report.{pop}.{chr}')
        with open(omega_plus_output_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line == '' or line.startswith('//'):
                    continue
                else:
                    try:
                        position, omega, start, end, valid = line.split()
                    except ValueError:
                        raise ValueError(f'Unexpected format in {omega_plus_output_file}\n{line}')
                    
                    position = int(float(position)) # convert to int
                    omega = float(omega)
                    start = int(float(start))-1 # convert to 0-based
                    end = int(float(end))
                    valid = bool(int(valid))

                    if valid:
                        # if (chr, start, end, omega) not in results:
                            # results.append((chr, start, end, omega))
                        results_pop.append((chr, start, end, omega, position)) # Collect all valid results
                        if omega >= cutoff_pop[pop]:  # Check against the cutoff
                            outliers.append((chr, start, end))
                    else:
                        if omega >= cutoff_pop[pop]:
                            print(f'Warning: invalid outlier at {chr}:{position}\t{start}-{end}')
        # Remove duplicates
        outliers = list(set(outliers))   # remove duplicates
        
        # Sort the outliers by start position within this chromosome
        outliers.sort(key=lambda x: x[1])    # Note: should sort chromosome first, sort -k1,1 -k2,2n
        for (chr, start, end) in outliers:
            output_file.write(f'{chr}\t{start}\t{end}\n')
    output_file.close()
    
    # Remove duplicates from results
    results_pop = list(set(results_pop))
    
    # Write all valid results
    with open(os.path.join(sub_output_dir, f'OmegaPlus.all.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        results_pop.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in results_pop]))
    
    # Sort the results by statistic
    results_pop.sort(key=lambda x: x[3], reverse=True)
    
    # Write the top 30% results
    top_30_results = results_pop[:int(0.3 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'OmegaPlus.top30per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_30_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_30_results]))
    
    # Write the top 20% results
    top_20_results = results_pop[:int(0.2 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'OmegaPlus.top20per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_20_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_20_results]))


    top_cut_off = [0.001, 0.01, 200, 500]
    for top in top_cut_off:
        if top < 1:
            tops = max(1,int(float(top * len(results_pop))))
        else:
            tops = int(top)

        # Write the top outliers
        tops_outliers = results_pop[:tops]
        # Sort the outliers by chromosome and by start position
        tops_outliers.sort(key=lambda x: (x[0], x[1]))
        with open(os.path.join(sub_output_dir, f'OmegaPlus.{top}outliers.{pop}.bed'), 'w') as f:
            for (chr, start, end, omega, position) in tops_outliers:
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
    merge_outliers('OmegaPlus', cutoff)    

    

# merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge.sh')
# with open(merge_script,'w') as f:
#     f.write(f"#!/bin/bash\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.out\n")
#     f.write(f"{load_bedtools}\n\n\n")
#     for pop in populations:
#         input_file = os.path.join(sub_output_dir, f'OmegaPlus.outliers.{pop}.bed')
#         output_file = os.path.join(sub_output_dir, f'OmegaPlus.outliers.{pop}.merged.bed')
#         gff_file = os.path.join(input_dir,gff_filename)
#         f.write(f"bedtools merge -i {input_file} > {output_file}\n")

#         # Get candidate genes
#         method = "OmegaPlus"
#         f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
# | awk '{{$2="{method}"; print}}' \\
# > {output_file}.genes.gff

# cat {output_file}.genes.gff | awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
# > {output_file}.genes.txt
# """)

# os.system(f"cd {sub_script_dir} && sbatch {merge_script}")
# print(wrap79(f"The following script has been submitted to slurm:"))
# print(wrap79(merge_script))



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))