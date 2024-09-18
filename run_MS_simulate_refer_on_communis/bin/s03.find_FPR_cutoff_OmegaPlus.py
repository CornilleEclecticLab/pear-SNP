#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.find_FPR_cutoff_OmegaPlus.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/09/16 11:30:42
# @Description:
#    

import datetime
import sys
import textwrap
import os
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


# Set the previous results directory
s02_folder = os.path.join(
    output_dir, 's02.generate_positive_selection_detect_OmegaPlus')


# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Read populations
populations = {}
with open(os.path.join(
        input_dir, 's01.sample_tab_population.txt'), 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        else:
            ind, pop = line.split()[0:2]
            populations[pop] = populations.get(pop, [])
            populations[pop].append(ind)


# Set replicates and sets numbers
replicates = 10
sets = 100


# Set FPR
FPR = 0.05


# Calculate the cutoff for each population
for pop in populations.keys():
    # Load the previous results
    pop_results_set_top = []
    pop_results = []
    for rep in range(1, replicates+1):
        file = os.path.join(
            s02_folder, f'OmegaPlus_Report.{pop}.REP{rep}')
        with open(file, 'r') as f:
            set_read = f.read().strip().split("//")
            for i in range(1, len(set_read)):
                set_results = []
                set_i = set_read[i].strip().split("\n")
                for line in set_i[1:]:
                    if line.startswith('#') or line == '':
                        continue
                    else:
                        omega = float(line.split()[1])
                        set_results.append(omega)
                        pop_results.append(omega)
                set_results.sort(reverse=True)
                pop_results_set_top.append(set_results[0])

    # Check the number of results
    if len(pop_results_set_top) != replicates*sets:
        raise ValueError(wrap79(
            f'Error: {pop} results number is not correct. You expected '
            f'{replicates*sets} but got {len(pop_results_set_top)}'))
    print(f'{pop} results number is {len(pop_results)}.')

    # Sort the results of top values of all sets
    pop_results_set_top.sort(reverse=True)
    cutoff_index = int(len(pop_results_set_top)*FPR) - 1
    cutoff = pop_results_set_top[cutoff_index]
    print(f'{pop} FPR 5% cutoff:=> {cutoff}\n')
    with open(os.path.join(sub_output_dir, f'{pop}.OmegaPlus.FPR005.cutoff.txt'), 'w') as f:
        f.write(f'{pop} Omega cutoff\n')
        f.write(f'{cutoff}\n')



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))