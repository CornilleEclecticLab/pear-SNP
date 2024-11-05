#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.calculate_Pi_Dxy_across_chromosome.py
# @Version  : 1.1.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/03/23 12:54:37
# @Description:
#     Calculate Pi and Fst across chromosome from the output of Pixy and read Dxy from Stacks 
#     results (If provided).
# @Update: 1.1.0 2024-04-16 11:36:27
#   1. NEW: Accept the group order in the input file.
#   2. NEW: Count N values in the output.
#   3. Change the output file path.

# @Update: 1.2.0 2024-04-29 17:46:33
#   1. NEW: Accept the pixy Fst results.

# @Update: 1.2.1 2024-05-10 
#   1. Finish editing in the last version.
#   2. Changed the output header.

# @Update: 2.0.0 2024-8-05
#   1. Update the function to calculate Fst from Pixy, now calculating the average Fst by SNP, not 
#      by window. Although the results are the same, the method is different.
#   2. Ignore ??

# @Update: 3.0.0 2024-10-24 15:25:05
#   1. IMPORTANT: fix a bug to count the number of differences and comparisons in the output file.
#   2. Output the total number of differences and comparisons in the output file.
#   3. Minor changes in the header of the output Pi by chr file.


import datetime
import argparse
import sys
import textwrap
import os
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "3.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output', script_basename)
sub_script_dir = os.path.join(bin_dir, script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters


def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Read the input file
parser = argparse.ArgumentParser(description=wrap79(f'''
    Calculate Pi and Fst across chromosome from the output of Pixy and 
    read Dxy from Stacks results (*.out, If provided).'''))
parser.add_argument('input',
                    nargs='+',
                    type=argparse.FileType('r'),
                    help=wrap79(f'''
                        The results files from Pixy to be read in, 
                        as well as the slurm *.out (summary) for running Stacks file.
                    '''))
parser.add_argument('-o', '--order',
                    help='The group order file.',
                    dest='order',
                    type=argparse.FileType('r'))

parser.add_argument('-s', '--sample_pop',
                    help=wrap79(f'''The sample and population name in the input file for s02. 
                    e.g. "../input/s02.sample_tab_population.txt"'''),
                    dest='sample_pop',
                    type=argparse.FileType('r'))

args = parser.parse_args()


# Read group order for the outputs
group_order = None
if args.order:
    group_order = [line.strip()
                   for line in args.order.read().strip().split('\n')]

# Read sample and population name for the outputs
pop_count = {}
if args.sample_pop:
    for line in args.sample_pop:
        sample, pop = line.strip().split()
        pop_count[pop] = pop_count.get(pop, 0) + 1


# Define the function to read diff and comparisons for Pi from Pixy
def read_Pi(dic, dic_chr, io):
    lines = io.readlines() 
    if lines[0].startswith('pop'):
        lines=lines[1:]
    for line in lines:
        ls = line.split()
        pop = ls[0]
        chr = ls[1]
        diff = int(ls[6] if ls[6] != 'NA' else 0)
        comparisons = int(ls[7] if ls[7] != 'NA' else 0)

        dic[pop] = dic.get(pop, [0, 0])
        dic[pop][0] += diff
        dic[pop][1] += comparisons

        dic_chr[pop] = dic_chr.get(pop, {})
        dic_chr[pop][chr] = dic_chr[pop].get(chr, [0, 0])
        dic_chr[pop][chr][0] += diff
        dic_chr[pop][chr][1] += comparisons

    return dic, dic_chr


# Define the function to read diff and comparisons for Dxy from Pixy
def read_Dxy(dic, dic_chr, io):
    lines = io.readlines()
    if lines[0].startswith('pop'):
        lines = lines[1:]
    for line in lines:
        ls = line.split()
        pop1 = ls[0]
        pop2 = ls[1]
        pops = (pop1, pop2)
        chr = ls[2]
        diff = int(ls[7] if ls[7] != 'NA' else 0)
        comparisons = int(ls[8] if ls[8] != 'NA' else 0)
        if diff == 0 and comparisons == 0:
            continue
        dic[pops] = dic.get(pops, [0, 0])
        dic[pops][0] += diff
        dic[pops][1] += comparisons
        dic_chr[pops] = dic_chr.get(pops, {})
        dic_chr[pops][chr] = dic_chr[pops].get(chr, [0, 0])
        dic_chr[pops][chr][0] += diff
        dic_chr[pops][chr][1] += comparisons
    return dic, dic_chr


# Define the function to read Fst from stacks
def read_Fst_stacks(dic, io):
    Fst_dict = dic
    records_start = False
    for line in io:
        line = line.strip()

        if line.startswith('Population pair divergence'):
            records_start = True
            continue

        if records_start and line == '':
            records_start = False

        if records_start:
            info = line.split(';')[0]
            print(info)
            pop1, pop2 = info.split(': mean Fst: ')[0].split('-')
            Fst = info.split(': mean Fst: ')[1]
            Fst_dict[(pop1, pop2)] = Fst
            Fst_dict[(pop2, pop1)] = Fst
    return Fst_dict


# Define the function to round the number
def round_num(num):
    return str(round(num, 6))  # Round to 6 decimal places


# Define the function to calculate Pi or Dxy
def cul_pi_dxy(diff, comparisons):
    if comparisons == 0:
        return 'NA'
    else:
        return round_num(diff / comparisons)


# Define the function to Read Fst from Pixy
def read_Fst_pixy(dic, dic_chr, io):
    lines = io.readlines()
    if lines[0].startswith('pop'):
        lines = lines[1:]
    for line in lines:
        ls = line.split()
        pop1 = ls[0]
        pop2 = ls[1]
        pops = (pop1, pop2)
        chr = ls[2]
        fst = float(ls[5]) if ls[5] != 'NA' else None
        no_snps = int(ls[6])

        dic[pops] = dic.get(pops, [])
        dic_chr[pops] = dic_chr.get(pops, {})
        dic_chr[pops][chr] = dic_chr[pops].get(chr, [])

        if fst is not None:
            dic[pops].append([fst, no_snps])
            dic_chr[pops][chr].append([fst, no_snps])

    return dic, dic_chr


# Define the function to calculate Fst from Pixy
def cul_Fst_pixy(list):
    if len(list) == 0:
        return 'NA'
    else:
        return round_num(sum([m[0]*m[1] for m in list]) / sum([n[1] for n in list]))


# Read the input files (results from Pixy and Stacks)
Pi_dic = {}
Pi_dic_chr = {}
Dxy_dic = {}
Dxy_dic_chr = {}
Fst_dic_Pixy = {}
Fst_dic_chr_Pixy = {}
Fst_dic_stacks = {}
for i in args.input:
    if i.name.endswith('_pi.txt'):  # Read Pi from Pixy
        Pi_dic, Pi_dic_chr = read_Pi(Pi_dic, Pi_dic_chr, i)

    elif i.name.endswith('_dxy.txt'):  # Read Dxy from Pixy
        Dxy_dic, Dxy_dic_chr = read_Dxy(Dxy_dic, Dxy_dic_chr, i)

    elif i.name.endswith('_fst.txt'):  # Read Fst from Pixy
        Fst_dic_Pixy, Fst_dic_chr_Pixy = read_Fst_pixy(Fst_dic_Pixy, Fst_dic_chr_Pixy, i)

    elif i.name.endswith('.out'):  # Read Fst from Stacks
        Fst_dic_stacks = read_Fst_stacks(Fst_dic_stacks, i)

    else:
        print('Unknown file type. Please check the file name.')
        print('The file name should end with "_pi.txt", "_dxy.txt", "_fst.txt" or ".out".')
        print('The unknown file name is:', i.name)
        sys.exit(1)


# Calculate Pi and Dxy and output the results into files
# Calculate Pi from Pixy
opi = open(os.path.join(output_dir, 'Pi.txt'), 'w')
opic = open(os.path.join(output_dir, 'Pi_chr.txt'), 'w')
opi.write('pop\tPi\tN\tsum_diffs\tsum_comparisons\n') if pop_count else opi.write('pop\tPi\n')
header = ''
group_order = group_order if group_order else list(Pi_dic.keys())
if len(group_order) != len(Pi_dic):
    raise ValueError(wrap79('The group order is not the same as'
                            'the populations in the input file.'))

for pop in group_order:
    N = pop_count.get(pop, 0)
    sum_diffs = Pi_dic[pop][0]
    sum_comparisons = Pi_dic[pop][1]
    Pi = cul_pi_dxy(sum_diffs, sum_comparisons)
    print(f'Pi for {pop} is {Pi}, N={N}')
    opi.write(f'{pop}\t{Pi}\t{N}\t{sum_diffs}\t{sum_comparisons}\n') if pop_count else opi.write(
        f'{pop}\t{Pi}\n')

    expected_header = header = 'pop\t' + \
        '\t'.join(Pi_dic_chr[pop].keys()) + '\n'
    if header == '':
        header = expected_header
        opic.write(header)
    elif header != expected_header:
        print('The chromosome is not the same in different populations.')
        sys.exit(1)

    Pi_chr = []
    for chr in list(Pi_dic_chr[pop].keys()):
        Pi = cul_pi_dxy(Pi_dic_chr[pop][chr][0], Pi_dic_chr[pop][chr][1])
        Pi_chr.append(Pi)
        print(f'Pi for {pop} in chromosome {chr} is {Pi}')

    opic.write(f'{pop}\t' + '\t'.join(Pi_chr) + '\n')
opi.close()
opic.close()


# Calculate Dxy from Pixy
odxy = open(os.path.join(output_dir, 'Dxy.txt'), 'w')
odxyc = open(os.path.join(output_dir, 'Dxy_chr.txt'), 'w')
odxy.write('\t' + '\t'.join(Pi_dic.keys()) + '\n')
Dxy_result = {}
for pop1 in group_order:
    content = pop1 + '\t'
    for pop2 in group_order:
        if pop1 == pop2:
            continue

        pops = (pop1, pop2) if (pop1, pop2) in Dxy_dic else (pop2, pop1)

        Dxy = cul_pi_dxy(Dxy_dic[pops][0], Dxy_dic[pops][1])
        print(f'Dxy between {pop1} and {pop2} is {Dxy}')
        content += Dxy + '\t'
        Dxy_result[(pop1, pop2)] = Dxy

        for chr in list(Dxy_dic_chr[pops].keys()):
            Dxy = cul_pi_dxy(Dxy_dic_chr[pops][chr]
                             [0], Dxy_dic_chr[pops][chr][1])
            print(
                f'Dxy between {pop1} and {pop2} in chromosome {chr} is {Dxy}')
            odxyc.write(f'{pop1}\t{pop2}\t{chr}\t{Dxy}\n')
    odxy.write(content + '\n')
odxy.close()
odxyc.close()


# Output the Dxy_pixy and Fst_stack mixed matrix
# Stacks
header = 'Pops\t'+'\t'.join(group_order) + '\n'
content = ''
for m in range(len(group_order)):
    pop1 = group_order[m]
    line = pop1
    for n in range(len(group_order)):
        pop2 = group_order[n]
        if m == n:    # pop1 == pop2, diagonal, skip
            line += '\t-'
            continue
        elif m > n:  # Lower triangle, for Dxy
            # line += '\t-'
            # continue
            pops = (pop1, pop2) if (pop1, pop2) in Dxy_dic else (pop2, pop1)
            Dxy = Dxy_result.get(pops, 'NA')
            line += '\t' + Dxy

        elif m < n:  # Upper triangle, for Fst
            pops = (pop1, pop2) if (
                pop1, pop2) in Fst_dic_stacks else (pop2, pop1)
            Fst = Fst_dic_stacks.get(pops, 'NA')
            line += '\t' + Fst
            print(f'{pop1}\t{pop2}\tFst_stacks\t{Fst}')
    line += '\n'
    content += line
open(os.path.join(output_dir, 'Fst_Stacks_in_Upper_and_Dxy_Pixy_in_Lower_matrix.tsv'),
     'w').write(header+content)


# Output the Dxy_pixy and Fst_pixy mixed matrix
# Pixy
header = 'Pops\t'+'\t'.join(group_order) + '\n'
content = ''
fst_chr_out = open(os.path.join(output_dir, 'Fst_Pixy_chr.txt'), 'w')
for m in range(len(group_order)):
    pop1 = group_order[m]
    line = pop1
    for n in range(len(group_order)):
        pop2 = group_order[n]
        if m == n:    # pop1 == pop2, diagonal, skip
            line += '\t-'
            continue
        elif m > n:  # Lower triangle, for Dxy
            pops = (pop1, pop2) if (pop1, pop2) in Dxy_dic else (pop2, pop1)
            Dxy = Dxy_result.get(pops, 'NA')
            line += '\t' + Dxy

        elif m < n:  # Upper triangle, for Fst
            pops = (pop1, pop2) if (
                pop1, pop2) in Fst_dic_Pixy else (pop2, pop1)
            Fst = cul_Fst_pixy(Fst_dic_Pixy.get(pops, []))
            line += '\t' + Fst
            print(f'{pop1}\t{pop2}\tFst_Pixy\t{Fst}')

            for chr in list(Fst_dic_chr_Pixy[pops].keys()):
                Fst_chr = cul_Fst_pixy(Fst_dic_chr_Pixy[pops].get(chr, []))
                print(f'{pop1}\t{pop2}\t{chr}\tFst_Pixy\t{Fst}')
    line += '\n'
    content += line
open(os.path.join(output_dir, 'Fst_Pixy_in_Upper_and_Dxy_Pixy_in_Lower_matrix.tsv'),
     'w').write(header+content)


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
