#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s02.choosing_pure_samples_manually.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/18 14:52:50
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
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
input_file_sample_list= os.path.join(work_dir,'output','s01.find_pure_individuals','s01.non_admix.id_map.txt')  
output_dir = os.path.join(work_dir,'output',script_basename)
output_list_for_bcftools = os.path.join(output_dir,'s02.non_admixed_list.txt')
output_nex_for_PAUP = os.path.join(output_dir,'s02.taxpartitions.nex')
output_txt_for_PAUP = os.path.join(output_dir,'s02.taxpartitions.txt')

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)

# Usage
print(wrap(f'''
Usage:  python3 {script_basename}.py 
'''))
print('')
print(wrap(f'''
Before running this script, you should have run s01.find_pure_individuals.py, and have manually edited the intput file {input_file_sample_list} ready.'''))
inputs = input('Please confirm that you have done the above steps. (y/n): ')
if inputs != 'y':
    print('Program terminated.')
    sys.exit(1)
print('')


nex_header = '''  BEGIN SETS;
  	TAXPARTITION SPECIES =
'''

nex_tail = '''  END;
'''


dict_species = {}
with open(input_file_sample_list,'r') as fi:
    for line in fi:
        if line.startswith('#'):
            continue
        else:
            line = line.strip()
            if line == '':
                continue
            else:
                line = line.split()
                id = line[0]
                sp = line[1]
                dict_species[sp] = dict_species.get(sp,[])
                dict_species[sp].append(id)


dict_num = {}
with open (output_txt_for_PAUP,'w') as fo:
    with open(output_list_for_bcftools,'w') as fo2:
        print(wrap(
            f'''Note:  The following individuals with index numbers are chosen as pure individuals, and will be useful for PAUP analysis setting up outgroup. The order of samples maybe CHANGED from the original id_map file <s01.non_admix.id_map.txt>. The content is saved in {output_txt_for_PAUP}:'''))

        n = 0
        for species,ids in dict_species.items():
            for id in ids:
                n += 1
                dict_num[id] = str(n)
                fo.write(f"{str(n)}\t{id}\t{species}\n")
                print(f"{str(n)}\t{id}\t{species}")
                fo2.write(id + '\n')


with open(output_nex_for_PAUP,'w') as fo:
    fo.write(nex_header)
    items = []
    for species in dict_species:
        n_min = dict_num[dict_species[species][0]]
        n_max = dict_num[dict_species[species][-1]]
        item = f"{species}: {n_min}-{n_max}"
        items.append(item)

    fo.write('  \t\t' +',\n  \t\t'.join(items) + ';\n')
    fo.write(nex_tail)


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))