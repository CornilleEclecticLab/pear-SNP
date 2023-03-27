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
from collections import defaultdict


start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))

message = \
'''The sequence ID in the alignment file that you assumed 
as root in the phylogenetic tree.'''

parser = argparse.ArgumentParser(
    prog="Remove_N",
    description="Remove N bases positions from root sequence in alignment")
parser.add_argument(
    'filename',
    metavar="filename.phy", help='phy Alignment file name')
parser.add_argument(
    '-r', '--rootID',
    help=message)

args = parser.parse_args()


# Read the alignment file.
with open(args.filename, 'r') as fo:
    seq_dic = defaultdict(str)
    first_line = True
    length = True
    max_len = 0

    for line in fo:
        if first_line:
            first_line = False
            continue
        line = line.strip().split()
        ID = line[0]
        seq = line[1]
        
        if len(ID) > max_len:
            max_len = len(ID)
        
        if length == True:
            length = len(seq)
        else:
            if len(seq) != length:
                print("ODD sequence length: " + ID)
        seq_dic[ID] = seq



# Calculate the rate of N on every site of assumed root sequences
remove_sites = []
if args.rootID is not None: 
    ids = args.rootID.strip().split(',')
    for i in range(length):
        count = sum(1 for id in ids if seq_dic[id][i] == 'N')
        rate = count / len(ids)
        if rate >= 0.6:
            remove_sites.append(i)


print("Input sequence length:{}".format(length))


# Remove the sites
# abridged_seq_dic = {id: [seq[i] for i in range(length) if i not in remove_sites] for id, seq in seq_dic.items()}


abridged_seq_dic = defaultdict(list)
for id, seq in seq_dic.items():
    for i in range(length):
        if i not in remove_sites:
            abridged_seq_dic[id].append(seq[i])

# abridged_seq_dic = {id: [seq[i] for i in range(length) if i not in remove_sites] for id, seq in seq_dic.items()}


# for id in seq_dic.keys():
#     abridged_seq_dic[id] = []

# for i in range(length):
#     if i not in remove_sites:
#         for id, seq in seq_dic.items():
#             abridged_seq_dic[id].append(seq[i])



# Generate the new file
with open(args.filename+".abridged.phy", 'w') as fo:
    fo.write("{} {}\n".format(len(seq_dic.keys()), length - len(remove_sites)))
    for id, seq in abridged_seq_dic.items():
        fo.write(id.ljust(max_len+4)+''.join(seq)+'\n')



end_time = datetime.datetime.now()
print('')
print(' END '.center(40, '='))
print(str(end_time-start_time).center(40))
