#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.catch_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024-09-18 16:43:01
# @Description:
#    Catch outliers in RAiSD output.

# Update: v2.0.0 2024-10-01 17:10:26
#    1. Get candidate genes from the outliers.

# Update: v2.1.0 2024-10-08 16:37:52
#    1. Use more cutoffs for outliers.
#    2. Write all valid results to a file.
# Update: v3.0.0 2025-01-23 00:09:28
#    1. Update output file name
#    2. Add 0.005 cutoff
#    3. Add 0.02 cutoff
# Update: v3.1.0 2025-02-13
#    1. Remove the duplicated genes from the output file.

import datetime
import sys
import textwrap
import os
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "3.1.0"
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


# Catch outliers from RAiSD by window/grid output, determined by the s00_config.py
# half_grid_window_size = grid_window_size // 2
pop_outliers_count = {}
for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'RAiSD.outliers.{pop}.bed'), 'w')
    outliers_count_pop = 0
    results_pop = []
    for chr in chromosomes:
        outliers = []
        # Read the RAiSD output
        raisd_output_file = os.path.join(output_dir,raisd_output_dir,raisd_files_name.format(pop=pop,chr=chr))
        with open(raisd_output_file) as f:
            # print(f'Reading {pop} {chr} {raisd_output_file}')
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
                    
                    # Get the window
                    if start == 0 and end == 0:
                        raise ValueError(
                            f'Unexpected start and end in {raisd_output_file}\n{line}')
                        # start = max(position - half_grid_window_size - 1, 0)
                        # end = position + half_grid_window_size
                    else:
                        start = start-1
                    
                    # Record the results
                    # Note: the following method to remove duplicates is very slow
                    # if (chr, start, end, mu) not in results:
                    #    results.append((chr, start, end, mu))
                    results_pop.append((chr, start, end, mu, position))
                    
                    # Check if it is an outlier
                    if mu >= cutoff_pop[pop]:
                        # Take care, this method to remove duplicates is very slow
                        # if (chr, start, end) not in outliers:
                            # outliers.append((chr, start, end))
                        outliers.append((chr, start, end))
                        # else:
                            # print(f'Outlier {chr}:{position} {start} {end} already found in {pop}')

        # Remove duplicates
        outliers = list(set(outliers))
        
        print(f'Outliers in {pop} {chr}: {len(outliers)}')

        # sort the outliers by start position within this chromosome
        if len(outliers) > 0:
            outliers.sort(key=lambda x: x[1])    # Note: should sort chromosome first, sort -k1,1 -k2,2n
            outliers_count_pop += len(outliers)
            for chr, start, end in outliers:
                output_file.write(f'{chr}\t{start}\t{end}\n')
    if outliers_count_pop == 0:
        print(f'No outlier was found in {pop}')
    output_file.close()
    
    # Remove duplicates from results
    results_pop = list(set(results_pop))
    
    # Write all valid results
    with open(os.path.join(sub_output_dir, f'RAiSD.all.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        results_pop.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{mu}' for chr, start, end, mu, position in results_pop]))
    
    # Sort the results by mu
    results_pop.sort(key=lambda x: x[3], reverse=True)
    
    # Write the top 30% results
    top_30_results = results_pop[:int(0.3 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'RAiSD.top30per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_30_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_30_results]))
    
    # Write the top 20% results
    top_20_results = results_pop[:int(0.2 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'RAiSD.top20per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_20_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chr}\t{position}\t{omega}' for (chr, start, end, omega, position) in top_20_results]))
    
    top_cut_off = [0.02, 0.01, 0.005, 0.001, 200, 500]
    for top in top_cut_off:
        if top < 1:
            tops = max(1,int(float(top * len(results_pop))))
        else:
            tops = int(top)
        tops_outliers = results_pop[:tops]
        tops_outliers.sort(key=lambda x: (x[0], x[1]))
        with open(os.path.join(sub_output_dir, f'RAiSD.{str(top).replace(".","_")}outliers.{pop}.bed'), 'w') as f:
            for chr, start, end, mu, position in tops_outliers:
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
            input_file = os.path.join(
                sub_output_dir, f'{method}.{str(cutoff).replace(".","_")}outliers.{pop}.bed')
            output_file = os.path.join(
                sub_output_dir, f'{method}.{str(cutoff).replace(".","_")}outliers.{pop}.merged.bed')
            gff_file = os.path.join(input_dir, gff_filename)
            f.write(f"bedtools merge -i {input_file} > {output_file}\n\n")

            # Get candidate genes
            f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
| awk '{{$2="{method}"; print}}' \\
> {output_file}.genes.gff

cat {output_file}.genes.gff \\
| awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
> {output_file}.genes.txt

if [[ -s "{output_file}.genes.txt" ]]; then
    sort -u {output_file}.genes.txt > {output_file}.genes.tmp && \\
    mv {output_file}.genes.tmp {output_file}.genes.txt
fi


""")

    os.system(f"cd {sub_script_dir} && sbatch {merge_script}")
    print(wrap79(f"The following script has been submitted to slurm:"))
    print(wrap79(merge_script))

ls = top_cut_off.copy()
ls.append('')
for cutoff in ls:
    merge_outliers('RAiSD', cutoff)


# # Merge the outliers
# merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge.sh')
# with open(merge_script,'w') as f:
#     f.write(f"#!/bin/bash\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
#     f.write(f"#SBATCH -o {os.path.basename(merge_script)}.%J.out\n\n")
#     f.write(f"{load_bedtools}\n\n\n")
#     for pop in populations:
#         input_file = os.path.join(sub_output_dir, f'RAiSD.outliers.{pop}.bed')
#         output_file = os.path.join(sub_output_dir, f'RAiSD.outliers.{pop}.merged.bed')
#         f.write(f"bedtools merge -d {grid_window_size} -i {input_file} > {output_file}\n\n")
#         # Get candidate genes
#         gff_file = os.path.join(input_dir, gff_filename)
#         method = "RAiSD"
#         f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
# | awk '{{$2="{method}"; print}}' \\
# > {output_file}.genes.gff

# cat {output_file}.genes.gff | awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
# > {output_file}.genes.txt
# """)

# os.system(f'cd {sub_script_dir} && sbatch {merge_script}')
# print(wrap79(f"The following script has been submitted to slurm:"))
# print(wrap79(merge_script))


# Catch outliers from RAiSD by window output
# Note: the following method is deprecated because not used cutoffs.
# threshold = 0.001  # Top 1% of the distribution

# for pop in populations:
#     # open output file
#     output_file = open(os.path.join(
#         sub_output_dir, f'RAiSD.outliers.by_win.{pop}.bed'), 'w')
#     outliers_pop = []
#     mu_count = 0
#     for chr in chromosomes:
#         results_chr = []
#         # Read the RAiSD output
#         raisd_output_file = os.path.join(
#             output_dir, raisd_output_dir, raisd_files_by_win_name.format(pop=pop, chr=chr))
#         with open(raisd_output_file) as f:
#             lines = f.readlines()
#             for line in lines:
#                 line = line.strip()
#                 if line == '':
#                     continue
#                 if line.startswith('//') or line.startswith('Position'):
#                     continue
#                 try:
#                     position, start, end, var, sfs, ld, mu = line.split()
#                 except ValueError:
#                     try:
#                         position, mu = line.split()
#                     except ValueError:
#                         raise ValueError(
#                             f'Unexpected format in {raisd_output_file}\n{line}')
                
#                 mu = float(mu)
#                 mu_count += 1
#                 start = int(float(start))
#                 end = int(float(end))
                
#                 if mu > 0:
#                     results_chr.append((chr,start, end, mu))
#                 else:
#                     raise ValueError(f'Unexpected mu value in {raisd_output_file}\n{line}')

#         if len(results_chr) > 0:
#             # Note: should sort chromosome first, sort -k1,1 -k2,2n
#             outliers_pop.extend(results_chr)

#     outliers_pop.sort(key=lambda x: x[3], reverse=True)
#     outliers_pop = outliers_pop[:int(len(outliers_pop)*threshold)]
#     # sort the outliers by start position
#     outliers_pop.sort(key=lambda x: (x[0], x[1]))

#     for chr, start, end, mu in outliers_pop:
#         output_file.write(f'{chr}\t{start}\t{end}\n')
#     output_file.close()

# # Merge the outliers
# merge_script = os.path.join(sub_script_dir, f'{script_basename}.merge_by_win.sh')
# with open(merge_script, 'w') as f:
#     f.write(f"#!/bin/bash\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.err\n")
#     f.write(f"#SBATCH -e {os.path.basename(merge_script)}.%J.out\n\n")
#     f.write(f"{load_bedtools}\n\n\n")
#     for pop in populations:
#         input_file = os.path.join(
#             sub_output_dir, f'RAiSD.outliers.by_win.{pop}.bed')
#         output_file = os.path.join(
#             sub_output_dir, f'RAiSD.outliers.by_win.{pop}.merged.bed')
#         f.write(
#             f"bedtools merge -d {grid_window_size} -i {input_file} > {output_file}\n\n")
#         # Get candidate genes
#         gff_file = os.path.join(input_dir, gff_filename)
#         method = "RAiSD"
#         f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
# | awk '{{$2="{method}"; print}}' \\
# > {output_file}.genes.gff

# cat {output_file}.genes.gff | awk -F'[;= ]' '{{for(i=1;i<=NF;i++) if($i=="Name") print $(i+1) "\t{method}"}}' \\
# > {output_file}.genes.txt
# """)

# os.system(f'cd {sub_script_dir} && sbatch {merge_script}')
# print(wrap79(f"The following script has been submitted to slurm:"))
# print(wrap79(merge_script))


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
