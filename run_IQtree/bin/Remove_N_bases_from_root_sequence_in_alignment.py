#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     :   Remove_N_bases_from_root_sequence_in_alignment.py
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/03/16 17:47:46
# Description:
#

import datetime
import argparse
import warnings

start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))

parser = argparse.ArgumentParser(prog="Remove_N",
                                 description="Remove_N_bases_from_root_sequence_in_alignment")
parser.add_argument('filename',
                    metavar="filename.phy", help='phy Alignment file name')
parser.add_argument('-r', '--rootID',
                    help='The sequence ID in the alignment file that you assumed as root in the phylogenetic tree')

args = parser.parse_args()


# Read the alignment file.
with open(args.filename, 'r') as fo:
    seq_dic = {}
    first_line = True
    length = True
    max = 0

    for line in fo:
        if first_line:
            first_line = False
            continue
        line = line.strip().split()
        ID = line[0]
        seq = line[1]
        
        if len(ID) > max:
            max = len(ID)
        
        if length == True:
            length = len(seq)
        else:
            if len(seq) != length:
                warnings.warn("ODD sequence length: " + ID)
        seq_dic[ID] = seq



# Calculate the rate of N on every site of assumed root sequences
remove_sites = []
ids = args.rootID.strip().split(',')
for i in range(length):
    count = 0

    for id in ids:
        if seq_dic[id][i] == 'N':
            count += 1

    rate = count / len(ids)
    if rate >= 0.6:
        remove_sites.append(i)

print(length)

# Remove the sites
abridged_seq_dic = {}
for id in seq_dic.keys():
    abridged_seq_dic[id] = []
    
for i in range(length):
    if i not in remove_sites:
        for id, seq in seq_dic.items():
            abridged_seq_dic[id].append(seq[i])


# Generate the new file
with open(args.filename+".abridged.phy", 'w') as fo:
    fo.write(str(len(seq_dic.keys()))+' '+str(length-len(remove_sites)) + '\n')
    for id, seq in abridged_seq_dic.items():
        fo.write(id.ljust(max+4)+''.join(seq)+'\n')


# map, zip, lambda用法


end_time = datetime.datetime.now()
print('')
print(' END '.center(40, '='))
print(str(end_time-start_time).center(40))
