#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.format_stacks_output.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/03/22 22:52:42
# @Description:
#     

# Update: v2.0.0 2024-10-08 10:23:57
#   1. Format the output of Stack populations to a table of statistics. 

import datetime
import sys
import os
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "2.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output',script_basename)


populations = []
with open(os.path.join(input_dir,'sample_tab_population.txt'),'r') as fo:
    for line in fo:
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        pop = line.split('\t')[1]
        if pop not in populations:
            populations.append(pop)


# Read data (expected Fst) as upper triangle matrix
Fst_dict = {}
with open(sys.argv[1],'r') as fi:
    records_start = False
    for line in fi:
        line = line.strip()

        if line.startswith('Population pair divergence'):
            records_start = True
            continue

        if records_start and line == '':
            records_start = False

        if records_start:
            info = line.split(';')[0]
            print(info)
            pop1,pop2 = info.split(': mean Fst: ')[0].split('-')
            Fst = info.split(': mean Fst: ')[1]
            Fst_dict[(pop1,pop2)] = Fst
            Fst_dict[(pop2,pop1)] = Fst


# Read data (Dxy) as upper triangle matrix



# Output the matrix
header = '\t'+'\t'.join(populations) + '\n'
content = ''
for m in range(len(populations)):
    pop1 = populations[m]
    line = pop1 + '\t'
    for n in range(len(populations)):
        pop2 = populations[n]
        if pop1 == pop2:
            line += '-\t'
            continue
        if n < m:
            line += '-\t'
            continue
        if n > m:
            if (pop1,pop2) in Fst_dict:
                line += Fst_dict[(pop1,pop2)] + '\t'
                print(f'{pop1}\t{pop2}\t{Fst_dict[(pop1,pop2)]}')
            else:
                line += 'NA\t'
                print(f'{pop1}\t{pop2}\tNA')
    line += '\n'
    content += line
open('The_Fst_Dxy_matrix.tsv','w').write(header+content)


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))