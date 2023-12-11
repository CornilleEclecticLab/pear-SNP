#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : choose_prefered_samples_in_a_group_of_clones.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/08/01 21:23:29
# @Description:
#

import datetime
import sys
import os
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



# Set path
version = "1.0.0"
script_basename = ''
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')

prefer_list = os.path.join(work_dir, 'input', 's04.prefered_sample_list.txt')
s04_plink2_out_dir = os.path.join(work_dir,'output', 'pear_Jul2023.variant_m2M2_maf005','pear_Jul2023.Combine_Chr.geno20_maf005.s04.clone_filtering')
in_list = os.path.join(s04_plink2_out_dir, 'without_clone0.354.king.cutoff.in.id')
out_list = os.path.join(s04_plink2_out_dir, 'without_clone0.354.king.cutoff.out.id')
kinship = os.path.join(s04_plink2_out_dir, 'without_clone0.354.kin0')
substitute_in_list = os.path.join(s04_plink2_out_dir, 'without_clone0.354.king.cutoff.substitute.in.id')
substitute_out_list = os.path.join(s04_plink2_out_dir, 'without_clone0.354.king.cutoff.substitute.out.id')


# Read in files
with open(prefer_list,'r') as f:
    prefer_list = f.read().strip().split('\n')

with open(in_list,'r') as f:
    in_list = f.read().strip().split('\n')[1:]
    in_list = [i.split('\t')[1] for i in in_list]

with open(out_list,'r') as f:
    out_list = f.read().strip().split('\n')[1:]
    out_list = [i.split('\t')[1] for i in out_list]


dic_kinship = {}
with open(kinship,'r') as f:
    f = f.read().strip().split('\n')[1:]
    for line in f:
        line = line.split('\t')
        ind1 = line[1]
        ind2 = line[3]
        kinship = float(line[7])

        if kinship >= 0.354:
            if ind1 in prefer_list:
                dic_kinship[ind1] = dic_kinship.get(ind1,[])
                dic_kinship[ind1].append(ind2)
            elif ind2 in prefer_list:
                dic_kinship[ind2] = dic_kinship.get(ind2,[])
                dic_kinship[ind2].append(ind1)



# Substitute the samples in out_list with the samples in prefer_list
no_substitution = list()
substitute_first = True
for i in range(len(out_list)):
    ind = out_list[i]
    if ind in prefer_list:
        substitute = False

        for j in dic_kinship[ind]:
            if j not in out_list and j not in prefer_list:

                if substitute_first:
                    print("Following samples are substituted with prefered samples:")

                substitute = True
                substitute_first = False

                print('\t', j , "\t->\t", out_list[i], sep='')
                in_list[in_list.index(j)] = out_list[i]
                out_list[i] = j
                break

        if not substitute:
            no_substitution.append(out_list[i])

if len(no_substitution) > 0:
    print('\nNo substitution for:\n\t', '\n\t'.join(no_substitution),sep='')



# Write the substituted in_list and out_list
with open(substitute_in_list,'w') as f:
    f.write('\n'.join(in_list))

with open(substitute_out_list,'w') as f:
    f.write('\n'.join(out_list))



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
