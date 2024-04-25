#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.random_chose_sample_and_set_haplotype_number.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/04/18 19:07:25
# @Description:
#     This script is used to randomly choose at most 9 individuals from each population and set haplotype number for each individual.
#     If the population has more than 9 individuals, randomly choose 9 individuals.
#     If the population has less than 9 but more than 3 individuals, choose all individuals.
#     This script would assign haplotype number for each individual, and set repeat number(s) for each individual.

# @Update: v1.1.0 2024-04-22 15:41:42
#     Generated new format output, that count haplotype number within each repeat. Because the Out of Memory error, I have to reduce the input vcf files by Repetition.

import datetime
import sys
import textwrap
import os
import random
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.1.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Read the population list
pop_dic = {}
with open(os.path.join(input_dir,'s01.individual_and_population_list.txt'),'r') as fi:
    for line in fi:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue
        ind, pop = line.split()
        pop_dic[pop] = pop_dic.get(pop, [])
        pop_dic[pop].append(ind)

# Randomly choose at most 9 individuals from each population
random_inds = {}
for pop, inds in pop_dic.items():
    if len(inds) >= 9:
        inds = random.sample(inds, 9)
    elif len(inds) >= 3:
        inds = random.sample(inds, len(inds))
        print(f'Warning: Population {pop} has {len(inds)} individuals, choose all individuals')
    elif len(inds) < 3:
        raise ValueError(f'No enough individuals in population:{pop}, '
                         f'you need at least 3 individuals in each population')
    elif len(pop) == 0:
        raise ValueError('No individual in population')

    random_inds[pop] = inds


# Assign haplotype number for each individual for each Repeat
rep_ind_pop_dic = {} # {REP1: {POP1: [ind1,ind2,ind3]}}
for pop, inds in random_inds.items():
    ind_n = 0
    for ind in inds:
        ind_n += 1
        if len(inds) == 9:
            rep_n = ind_n % 3 if ind_n % 3 != 0 else 3
        else:
            times = 9 // len(inds)
            remain = 9 % len(inds)
            use_times = times + 1 if ind_n <= remain else times

            if use_times == 1:
                rep_n = ind_n % 3 if ind_n % 3 != 0 else 3

            elif use_times == 2:
                rep_n = str(ind_n % 3 if ind_n % 3 != 0 else 3) + ',' + str(ind_n % 3 + 1)

            elif use_times == 3:
                rep_n = '1,2,3'

        rep_n = str(rep_n)
        reps = [str(i) for i in rep_n.split(',')]
        for rep in reps:
            rep = 'REP' + rep
            rep_ind_pop_dic[rep] = rep_ind_pop_dic.get(rep, {})
            rep_ind_pop_dic[rep][pop] = rep_ind_pop_dic[rep].get(pop, [])
            rep_ind_pop_dic[rep][pop].append(ind)

# Write the random individuals to file
output_file = os.path.join(input_dir,'s01.random_individual_population.txt')
if os.path.exists(output_file):
    print("Output file already exists. Skipping writing operation.")
    exit(0)
else:
    with open(output_file,'w') as fo:
        fo.write("#Indvidual\tPopulation\tIndvidual_number\tHaplotype1_number\tHaplotype2_number\tRepeat_number\n")
        for rep, pop_inds in rep_ind_pop_dic.items():
            hap_n = 0
            ind_n = 0
            for pop, inds in pop_inds.items():
                for ind in inds:
                    hap1_n = hap_n
                    hap2_n = hap_n + 1
                    hap_n += 2
                    ind_n += 1
                    fo.write(f"{ind}\t{pop}\t{ind_n}\t{hap1_n}\t{hap2_n}\t{rep}\n")







# # Write the random individuals to file
# output_file = os.path.join(input_dir,'s01.random_individual_population.txt')
# if os.path.exists(output_file):
#     print("Output file already exists. Skipping writing operation.")
# else:
#     with open(output_file),'w') as fo:
#         fo.write("#Indvidual\tPopulation\tIndvidual_number\tHaplotype1_number\tHaplotype2_number\tRepeat_number\n")
#         hap_n = 0
#         for pop, inds in random_inds.items():
#             ind_n = 0
#             # check_rep = {}
#             for ind in inds:
#                 hap1_n = hap_n 
#                 hap2_n = hap_n + 1
#                 hap_n += 2
#                 ind_n += 1
#                 if len(inds) == 9:
#                     rep_n = ind_n % 3 if ind_n % 3 != 0 else 3
#                     # check_rep[rep_n] = check_rep.get(rep_n, 0) + 1

#                 else:
#                     times = 9 // len(inds)
#                     remain= 9 % len(inds) 
#                     use_times = times + 1 if ind_n <= remain else times

#                     if use_times == 1:
#                         rep_n = ind_n%3 if ind_n%3 != 0 else 3
#                         # check_rep[rep_n] = check_rep.get(rep_n, 0) + 1

#                     elif use_times == 2:
#                         rep_n = str(ind_n%3 if ind_n%3!=0 else 3)+','+str(ind_n%3+1)
#                         # check_rep[rep_n] = check_rep.get(rep_n, 0) + 1

#                     elif use_times == 3:
#                         rep_n = '1,2,3'
#                         # check_rep[rep_n] = check_rep.get(rep_n, 0) + 1

#                 rep_n = str(rep_n)
#                 rep = 'REP' + ',REP'.join([str(i) for i in rep_n.split(',')])
#                 fo.write(f"{ind}\t{pop}\t{ind_n}\t{hap1_n}\t{hap2_n}\t{rep}\n")

#             # # For Developing: Check the repeat number 
#             # for rep, count in check_rep.items():
#             #     if count > 3:
#             #         raise ValueError(f'Repeat {rep} has more than 3 individuals in population {pop}')
#             #     for r in str(rep).split(','):
#             #         if int(r) > 3:
#             #             raise ValueError(f'Population {pop} has more than 3 repeats,{rep}')




end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
