#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.combine_cross_coal.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/04/19 19:52:59
# @Description:
#
# @History:
#    v1.0.1 2024-04-29 14:20:37 bug fix

import datetime
import sys
import textwrap
import os
import glob
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.0.1"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output', script_basename)
sub_script_dir = os.path.join(bin_dir, script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)
os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters


def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Define the msmc2 output directory
msmc2_within_pop_dir = os.path.join(
    work_dir, 'output', 's04.run_MSMC2_within_population')
msmc2_across_pop_dir = os.path.join(
    work_dir, 'output', 's05.run_MSMC2_across_population')


# Function to get path to the msmc2 output file
def get_msmc2_out_path(pop1, rep1, *args):
    if len(args) == 0:
        return os.path.join(msmc2_within_pop_dir, f'{pop1}.{rep1}.msmc2.final.txt')
    elif len(args) == 2:
        pop2, rep2 = args
        return os.path.join(msmc2_across_pop_dir, f'{pop1}.{rep1}.{pop2}.{rep2}.msmc2.final.txt')
    else:
        raise ValueError('Invalid number of arguments')


# Write the script to combine the msmc2 output files
with open(os.path.join(sub_script_dir, script_basename+'.sh'), 'w') as fo:
    fo.write(f"""#!/usr/bin/env bash
# This script is generated ../{script_basename}.py by {version} on {start_time}
#SBATCH -J {script_basename}
#SBATCH -o {script_basename}.%J.out
#SBATCH -e {script_basename}.%J.err
#SBATCH --mem=1G

source {os.path.join(bin_dir,'s00.load_environment.sh')}
""")
# Glob the msmc2 output files (based on the across population output)
    count_msmc2_files = 0
    for path in glob.glob(os.path.join(msmc2_across_pop_dir, '*.msmc2.final.txt')):
        basename = os.path.basename(path)
        pop1, rep1, pop2, rep2 = basename.split('.')[0:4]
        count_msmc2_files += 1
        print(f'{count_msmc2_files} Processing {pop1} {rep1} {pop2} {rep2} {basename}')
        basename_no_ext = os.path.splitext(basename)[0]
        fo.write(f"""
combineCrossCoal.py \\
    {get_msmc2_out_path(pop1,rep1,pop2,rep2)} \\
    {get_msmc2_out_path(pop1,rep1)}  \\
    {get_msmc2_out_path(pop2,rep2)}  \\
> {os.path.join(output_dir, basename_no_ext + '.combine.txt')}
""")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
