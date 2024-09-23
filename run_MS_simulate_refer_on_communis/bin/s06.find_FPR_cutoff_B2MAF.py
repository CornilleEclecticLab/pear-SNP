#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.find_FPR_cutoff_OmegaPlus.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024-09-23 16:12:11
# @Description:
#    

import datetime
import sys
import textwrap
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


# Set FPR
FPR = 0.05


# Set the previous results directory
s05_folder = os.path.join(
    output_dir, 's05.generate_balancing_selection_detect_BalLeRMixP')


# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)



# Calculate the cutoff for each population
results = glob.glob(os.path.join(s05_folder, '*.B2MAF.txt'))
population_results = {}
for result in results:
    pop = os.path.basename(result).split('.')[2]
    population_results[pop] = population_results.get(pop, [])
    population_results[pop].append(result)

for pop, results in population_results.items():
    pop_results = []
    print(f'Processing {pop}...')
    c=0
    s=0
    for result in results:
        clrs = []
        with open(result, 'r') as f:
            for line in f:
                if line.startswith('#') or line == '':
                    continue
                elif line.startswith('physPos'):
                    continue

                physPos, genPos, clr, x_hat, s_hat, A_hat, nSites = \
                    line.strip().split()
                s_hat = float(s_hat)
                s+=1
                if s_hat <= 1:
                    # print(line)
                    c+=1
                    continue
                elif s_hat > 1:
                    clr = float(clr)
                    clrs.append(clr)
        clrs.sort(reverse=True)
        best_clr = clrs[0]
        pop_results.append(best_clr)
    print(c,s)
    pop_results.sort(reverse=True)
    cutoff = pop_results[int(float(len(pop_results)*FPR)) - 1]
    print(f'cutoff for {pop} is {cutoff}\n')
    with open(os.path.join(
        sub_output_dir, f'{pop}.B2MAF.FPR005.cutoff.txt'), 'w') as f:
        f.write(f'{pop}\tB2MAF\tcutoff\n\
                  {cutoff}')




end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))