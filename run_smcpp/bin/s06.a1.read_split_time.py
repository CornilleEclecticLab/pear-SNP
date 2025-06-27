#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.a1.read_split_time.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024-08-16 18:10:01
# @Description:
#    This script is extract split time form the smc++ csv files.


import datetime
import os
import sys
import textwrap
from warnings import warn
import argparse
import csv
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


default_spline = ("piecewise", "cubic", "pchip")[0]  # piecewise, cubic, pchip
parser = argparse.ArgumentParser(
    description='Generate sub-scripts for the "smc++ estimate" command.')
parser.add_argument('-s', '--spline',
                    type=str,
                    default=default_spline,
                    choices=['piecewise', 'cubic', 'pchip'],
                    help='The type of spline to use for the analysis, "piecewise","cubic", or "pchip", default: %(default)s')
args = parser.parse_args()
spline_type = args.spline
print(f'The spline type is: {spline_type}')


version = "1.0.0"

s06_script_basename = "s06.plot_split_time_joint"+f'_{spline_type}'

script_basename = 's06.a1.read_split_time'+f'_{spline_type}'

script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = work_dir+'/input/'
output_dir = work_dir+'/output/'+script_basename
sub_script_dir = work_dir+'/bin/'+script_basename
s06_output_dir = work_dir+'/output/'+s06_script_basename
pop_order_file = input_dir+'s06.population_order.txt'

# Format: individual_name population_name
individual_population_list = input_dir + \
    "/s01.individuals_and_populations_list.txt"
# Format: population1 population2
population_pair_list = input_dir+"/s04.population_pair_list.txt"
load_singularity = '# module load system/singularity-3.7.3 # singularity is installed and callable on the cluster without loading a module'

cpu_cores = "20"
em_iterations = "50"


os.system("mkdir -p "+output_dir)


def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


population_pair_list = input_dir+"/s04.population_pair_list.txt"
pop_list = []
with open(population_pair_list, 'r') as f:
    if f.readline().strip()[0].startswith('#'):
        population_pairs = f.readlines()
    else:
        f.seek(0)
        population_pairs = f.readlines()

    population_pairs = [i.strip().split() for i in population_pairs]

    for pops in population_pairs:
        for pop in pops:
            pop_list.append(pop) if pop not in pop_list else None


if not os.path.exists(pop_order_file):
    warn(f"Warning: The file {pop_order_file} does not exist."
         f"Using the order in the population_pair_list.")
else:
    with open(pop_order_file, 'r') as f:
        pop_list = f.read().strip().split('\n')

dic_split_time = {}
for pops in population_pairs:

    for i in range(2):
        pop1, pop2 = pops[i], pops[i-1]
        print(f'{" "+pop2+"<<<from<<<"+pop1+" ":=^79}')
        with open(s06_output_dir+'/'+s06_script_basename+f".{pop1}.{pop2}.csv", 'r') as f: # read the csv file
            csv_reader = csv.DictReader(f)
            dic = {}
            for row in csv_reader:
                pop = row['label']
                time = float(row['x'])
                time = int(time)
                dic[pop] = dic.get(pop, 0)
                if time > dic[pop]:
                    dic[pop] = time

            if dic[pop1] > dic[pop2]:
                pass
            else:
                raise ValueError(f"Error: {pop1} time is earlier than {pop2} split time.")

            dic_split_time[(pop1, pop2)] = dic[pop2]
            print(f'{" "+str(dic[pop2])+" ": ^79}')


content = ''
for m in range(-1,len(pop_list)):
    for n in range(-1, len(pop_list)):
        pop_m = pop_list[m]
        pop_n = pop_list[n]

        if m == n:
            if m ==n == -1:
                content += 'From\To'
            else:
                content += '\t-'
        elif m == -1:
            if n > -1:
                content += '\t'+pop_n
        elif n == -1:
            if m > -1:
                content +=  pop_m
        elif m<n: # upper triangle
            #content += "\t\t+"
            content += '\t'+str(dic_split_time[(pop_m, pop_n)])
        elif m>n: # lower triangle
            #content += "\t\t+"
            content += '\t'+str(dic_split_time[(pop_m, pop_n)])
    content += '\n'

with open(output_dir+'/'+script_basename+'.split_time.tsv', 'w') as f:
    f.write(content)
