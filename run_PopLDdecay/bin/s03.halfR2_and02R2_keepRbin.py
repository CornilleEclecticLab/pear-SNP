#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.halfR2_and02R2_keepRbin.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/10/29 14:12:04
# @Description:
#    

import glob
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
input_dir = os.path.join(output_dir, "s01.generate_PopLDdecay_running_scripts")
sub_output_dir = os.path.join(output_dir,script_basename)
output_file_path = os.path.join(sub_output_dir, "results.txt")

# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Get list of .bin files in the input directory
file_list = glob.glob(f"{input_dir}/*.bin")
print(wrap79(
    f"Found {len(file_list)} PopLDdecay output files in {input_dir}\n"))


# Open output file for writing and write file list
ou = open(output_file_path, "w")
ou.write("PopLDdecay output file list:\n")
ou.write('\n'.join(file_list)+'\n\n')


# Get the distance when R^2 value is half of the maximum R^2 value
ou.write("pop\tdistance_bp_when_r2=1/2maxr2\n")
for file in file_list:
    with open(file, "r") as fh_bgz_in:
        # Initialize the first Dist and R^2 values
        dist_1 , r_sq1 = False, False

        for line in fh_bgz_in:
            # Skip the first line
            if line.startswith("#"):
                if line.strip().split()[1] != "Mean_r^2":
                    raise ValueError(
                        f"Expected 'Mean_r^2' in the second column of the first line of PopLDdecay result, "
                        f"but got '{line.strip().split()[1]}'"
                    )
                continue
                
            line = line.strip()
                
            if (dist_1 is False) and (r_sq1 is False):
                dist_1, r_sq1 = line.split("\t")[0:2]

            dist, r_sq = line.split("\t")[0:2]
                
            if (float(r_sq) <= (float(r_sq1) / 2)):
                pop = '_'.join(os.path.basename(
                    file).split(".")[0].split("_")[2:])
                ou.write(f"{pop}\t{dist}\n")
                break  # Exit loop after the first match
ou.write("\n")


# Get the distance when R^2 value is less than 0.2
ou.write("pop\tdistance_bp_when_r2<0.2\n")
for file in file_list:
    with open(file, "r") as fh_bgz_in:
        for line in fh_bgz_in:
            # Skip the first line
            if line.startswith("#"):
                continue
            # Initialize the first Dist and R^2 values
            dist_1, r_sq1 = False, False
        
            line = line.strip()
            
            if (dist_1 is False) and (r_sq1 is False):
                dist_1, r_sq1 = line.split("\t")[0:2]    

            dist, r_sq = line.split("\t")[0:2]
            if float(r_sq) < 0.2:
                pop = '_'.join(os.path.basename(
                    file).split(".")[0].split("_")[2:])
                ou.write(f"{pop}\t{dist}\n")
                break  # Exit loop after the first match


print(wrap79(f"Output written to {output_file_path}"))


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
