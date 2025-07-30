#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.catch_outliers.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025-07-21 10:34:37
# @Description:
#    Catch the outliers from pixy TajimaD output.

# Update: v1.0.0 
#   1. Modified from s04.catch_outliers.py for OmegaPlus v3.1.0
#   2. Ouput cutoff statistic values.


import datetime
import sys
import textwrap
import os
import numpy as np
from s00_config import *
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0" # For pyrifolia
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

method = 'Pixy_TajimaD'

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
# cutoff_pop = {}
# for pop in populations:
#     cutoff_file = os.path.join(input_dir,cutoff_files_dir,f'{pop}.{cutoff_files_suffix}')
#     with open(cutoff_file) as f:
#         lines = f.readlines()
#         header = lines[0].strip()
#         cutoff = lines[-1].strip()
#         if header.startswith('tajima_ds') or header.startswith('#') or header.startswith(pop):
#             cutoff_pop[pop] = float(cutoff)
#         else:
#             raise ValueError(f'Unexpected format in {cutoff_file}')


# Read TajimaD results from Pixy, and catch the outliers based on neutral envlope
pop_outliers_count = {}
pop_results_dict = {}
for pop in populations:
    # open output file
    output_file = open(os.path.join(
        sub_output_dir, f'{method}.outliers.{pop}.bed'), 'w')
    results_pop = []
    for chromosome in chromosomes:
        outliers = []
        # Read the Pixy_TajimaD output
        pixy_output_file = os.path.join(output_dir,pixy_output_dir,f'Pixy.{pop}.{chromosome}.1k_tajima_d.txt')

        f = open(pixy_output_file, 'r')
        lines = f.readlines()
        if lines[0].startswith('pop'):
            header = lines[0].strip()
        else:
            raise ValueError(f'Unexpected format in {pixy_output_file}\n{lines[0]}')
            
        for line in lines[1:]:
            line = line.strip()
            if line == '':
                continue

            try:
                pop, chromosome, window_pos_1, window_pos_2, tajima_d, no_sites, raw_pi, raw_watterson_theta, tajima_d_stdev = line.split()

            except ValueError:
                raise ValueError(f'Unexpected format in {pixy_output_file}\n{line}')
                    

            tajima_d = float(tajima_d)
            start = int(window_pos_1)-1 # convert to 0-based
            end = int(window_pos_2)
            position = int(0.5*(start+end)) # convert to int
                
            if tajima_d == 'NA':
                continue

            results_pop.append((chromosome, start, end, tajima_d, position)) # Collect all valid results
            # if tajima_d <= cutoff_pop[pop]:  # Check against the cutoff
            #     outliers.append((chromosome, start, end))

        # Remove duplicates
        outliers = list(set(outliers))   # remove duplicates
        
        # Sort the outliers by start position within this chromosome
        outliers.sort(key=lambda x: x[1])    # Note: should sort chromosome first, sort -k1,1 -k2,2n
        for (chromosome, start, end) in outliers:
            output_file.write(f'{chromosome}\t{start}\t{end}\n')
    output_file.close()
    
    # Remove duplicates from results
    results_pop = list(set(results_pop))
    
    # Write all valid results
    with open(os.path.join(sub_output_dir, f'{method}.all.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        results_pop.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chromosome}\t{position}\t{tajima_d}' for (chromosome, start, end, tajima_d, position) in results_pop]))
    
    # Sort the results by statistic
    # Sort by TajimaD value, negative TajimaD values as outliers
    results_pop.sort(key=lambda x: x[3], reverse=False)
    
    # Write the top 20% results
    top_20_results = results_pop[:int(0.2 * len(results_pop))]
    with open(os.path.join(sub_output_dir, f'{method}.top20per.{pop}.txt'), 'w') as f:
        # Sort by chromosome and position
        top_20_results.sort(key=lambda x: (x[0], x[1]))
        f.write('\n'.join([f'{chromosome}\t{position}\t{tajima_d}' for (chromosome, start, end, tajima_d, position) in top_20_results]))


    # top_cut_off = [0.01, 0.005, 0.001, 0.0001, 200, 500]
    top_cut_off = [f"{i:.4f}" for i in np.arange(0.05, 0.0009, -0.0005)]
    for top in top_cut_off:
        if float(top) < 1:
            tops = max(1,int(float(float(top) * len(results_pop))))
        else:
            tops = int(top)

        # Wirte the cutoff value
        cutoff_stat_value = results_pop[tops-1][3]

        with open(os.path.join(sub_output_dir, f'{pop}.{method}.{str(top).replace(".","_")}.cutoff.txt'), 'w') as f_cutoff:
            f_cutoff.write(f'{pop} {method} cutoff\n')
            f_cutoff.write(f'{cutoff_stat_value}\n')

        # Write the top outliers
        tops_outliers = results_pop[:tops]

        # Sort the outliers by chromosome and by start position
        tops_outliers.sort(key=lambda x: (x[0], x[1]))
        with open(os.path.join(sub_output_dir, f'{method}.{str(top).replace(".","_")}outliers.{pop}.bed'), 'w') as f:
            for (chromosome, start, end, d, position) in tops_outliers:
                f.write(f'{chromosome}\t{start}\t{end}\n')
    
    tajima_d_values_pop = np.array([d for chr, start, end, d, position in results_pop])
    zscores_pop = (tajima_d_values_pop - np.mean(tajima_d_values_pop)) / np.std(tajima_d_values_pop)
    results_with_zscores_pop = [
        (*item, z) for item, z in zip(results_pop, zscores_pop)
    ]

    
    # Add the results to the all pop results dict
    pop_results_dict[pop] = results_pop
    
    # Set a lot of cutoffs for Z-scores
    zscores_cutoff_list = [f"{i:.1f}" for i in np.arange(2.0, 10.01, 0.1) ]  # 2.0, 2.1, ... 10.0
    print(f"Z-scores cutoffs: {zscores_cutoff_list}")
    for z_cutoff in zscores_cutoff_list:
        for i, item in enumerate(results_with_zscores_pop):
            if abs(item[5]) < float(z_cutoff):
                ztops_outliers = results_with_zscores_pop[:i]
                break

        ztops_outliers.sort(key=lambda x: (x[0], x[1]))

        with open(os.path.join(sub_output_dir, f'{method}.Z{str(z_cutoff).replace(".","_")}outliers.{pop}.bed'), 'w') as f:
            for chr, start, end, mu, position, z in ztops_outliers:
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
        f.write(f"source {bin_dir}/s00.load_bedtools.sh\n\n\n")
        for pop in populations:
            input_file = os.path.join(
                sub_output_dir, f'{method}.{str(cutoff).replace(".","_")}outliers.{pop}.bed')
            output_file = os.path.join(
                sub_output_dir, f'{method}.{str(cutoff).replace(".","_")}outliers.{pop}.merged.bed')
            gff_file = os.path.join(input_dir, gff_filename)
            f.write(f"bedtools merge -i {input_file} > {output_file}\n\n")

            # Get candidate genes
            f.write(f"""bedtools intersect -a {gff_file} -b {output_file} -wa \\
| awk -v OFS='\\t' '{{$2="{method}"; print}}' \\
> {output_file}.genes.gff

cat {output_file}.genes.gff \\
| awk -F'[\\t;= ]' '{{for(i=1;i<=NF;i++) if($i=="Accession") print $(i+1) "\t{method}"}}' \\
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
# ls.append('') # This is from the cutoff control file
ls.extend([f'Z{z_cutoff}' for z_cutoff in zscores_cutoff_list])
for cutoff in ls:
    merge_outliers(method, cutoff)    


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))