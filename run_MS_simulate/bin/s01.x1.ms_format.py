#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.x1.ms_format.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/01/30 10:26:52
# @Description:
#    
# Format the ms like output file to print the formated output by: 
#  1. removing the duplicate positions;
#  2. adding a minor increment to the positions to avoid the float precision error in OmegaPlus and SweeD.


import sys
import os


version = "1.0.0"

# Usage
usage = f'''
    Usage:
    {sys.argv[0]} "ms_file" > "formatted_ms_file"
'''

# Print usage if needed
if len(sys.argv) != 2:
    print(usage)
    sys.exit(1)

# Control to output the duplicate positions or formated ms file
switch_count_duplicate_position = False
if switch_count_duplicate_position:
    duplicate_position_count = []


# Read ms/scrm output file and print the formatted file
scrm_file = sys.argv[1]

# Check the file exists
if not os.path.exists(scrm_file):
    print(f'Error: {scrm_file} not found')
    sys.exit(1)


# Function to check the positions
def check_float_precision_error(last_num, current_num, multipler):    
    int1 = int(last_num * multipler)
    int2 = int(current_num * multipler)

    if int1 < int2:  # pass the check
        return True
    elif int1 == int2:  # check the decimal part
        return False
    elif int1 > int2:
        raise ValueError(
            f'Error: the No {segment_count} positions not sorted')


# Read file
file = open(scrm_file, 'r')
segments = file.read().split("//")
segment_count = 0


# Get the accuracy of the positions
pars = segments[0].strip().splitlines()
for line in pars:
    if line.startswith("scrm") or line.startswith("ms"):
        parameters = line.split()
        read_parameters = True

        accuracy = int(parameters[parameters.index('-p') + 1])
        increment = 10 ** -accuracy
        minor_increment = increment / 10
        higher_accuracy = accuracy + 1

        nsites = int(parameters[parameters.index('-r') + 2])
        while minor_increment*nsites > 0.11:  # Ensure not shift the positions to the next integer
            minor_increment /= 10
            higher_accuracy += 1

        break
print(segments[0].strip())

# Read the segments
for seg in segments[1:]:
    lines = seg.strip().splitlines()
    # Count the number of segments
    segment_count += 1
    if switch_count_duplicate_position:
        duplicate_position_count_per_seg = 0

    for n, line in enumerate(lines):
        line = line.strip()
        # Get the number of segsites
        if line.startswith("segsites:"):
            segsites = int(line.split()[1])
            continue

        if line.startswith("positions:"):
            positions = [float(pos) for pos in line.split()[1:]]
            unique_positions = []
            unique_indices = []
            unique_variants = []
            for i, pos in enumerate(positions):
                pos += minor_increment  # Add a minor increment to avoid the float precision error in OmegaPlus and SweeD
                if i == 0:
                    unique_positions.append(pos)
                    unique_indices.append(i)
                    continue
                
                last_unique_pos = unique_positions[-1]
                
                if last_unique_pos > pos:
                    print(f"last_unique_pos: {last_unique_pos}")
                    raise ValueError(
                        f'Error: the No. {segment_count} positions not sorted')
                elif last_unique_pos == pos:
                    if switch_count_duplicate_position:
                        duplicate_position_count_per_seg += 1
                    continue
                elif last_unique_pos < pos:
                    # unique_positions.append(pos)
                    # unique_indices.append(i)
                    if check_float_precision_error(last_unique_pos, pos, nsites):
                        unique_positions.append(pos)
                        unique_indices.append(i)
                    else:
                        if switch_count_duplicate_position:
                            duplicate_position_count_per_seg += 1
                        raise ValueError(
                            f'Error: Not pass the check in the segment {segment_count}, with positions {last_unique_pos} and {pos}')
            if switch_count_duplicate_position:
                print(
                    f"Found {duplicate_position_count_per_seg} duplicate positions of {segsites} segsites in segment {segment_count}.")

            # Print the segsites and positions
            if not switch_count_duplicate_position:
                print('\n//')
                unique_segsites = len(unique_positions)
                print(f"segsites: {unique_segsites}")
                print(
                    f'positions: {" ".join(f"{float(pos):.{higher_accuracy}f}" for pos in unique_positions)}')
    
        # Deal with the haplotypes
        elif line.startswith("1") or line.startswith("0"):
            if not switch_count_duplicate_position:
                haplotype = [line[u] for u in unique_indices]
                print(f"{''.join(haplotype)}")
            pass
        else:
            raise ValueError(
                f'Error: Unknown line {n} in the segment {segment_count}')
file.close()
