#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s01.find_pure_individuals.py
# @Version  : 3.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/09/15 11:35:21
# @Description:
#     

# Update v2.0.0: 2023-09-20 12:48:24
# 1. Fix a bug that index of ids is not correct

# Update v2.0.1: 2023-09-21 13:28:33
# 1. Output the Q file with header and id_map
# 2. Output the best cluster for each individual in the id_map.txt file for further manual check

# Update v2.0.2: 2023-09-22 11:11:28
# 1. Remain the uni_id, and add species name column in the output file
# 2. Rename some variables

# Update v3.0.0: 2024-02-08 11:12:26


import argparse
import datetime
import os
import sys
import textwrap
from warnings import warn


start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


# Constants
version = "3.0.0"
threshold = 0.8 #0.8   # The threshold value that to distinguish the admixture, which is the maximum value of Qs for each individual to be considered as non-admixture. This value should be set according to the Q file.


# File paths
# script_basename = 's01.find_pure_individuals'
# script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
# work_dir = os.path.dirname(script_path)
# input_dir = os.path.join(work_dir,'input')
# input_file_Q = os.path.join(input_dir,'test.Q')                #######################################################
# input_file_fam = os.path.join(input_dir,'test.fam')            # Please change these paths to your actual input file, which is an output file from Plink and contains the samples ID #
# input_file_id_map = os.path.join(input_dir,'test.id_map')      # The id_map file is the sample ID in vcf file mapping onto the short ID (indicate species name and unique number) #
# output_dir = os.path.join(work_dir,'output',script_basename)
# output_file = os.path.join(output_dir,'s01.non_admix.id_map.txt')
# output_Q = os.path.join(output_dir,'s01.id_map_header.Q.tsv')
# sub_script_dir = os.path.join(work_dir,'bin',script_basename)

# os.system(f'mkdir -p {output_dir}')

def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--color_pallet')
parser.add_argument('-i', '--id_map', required=True)
parser.add_argument('-o', '--order_group')
parser.add_argument('-g', '--group_map')
parser.add_argument('-q', '--Q_file', action='extend' ,nargs='+')
parser.add_argument('-f', '--fam', action='extend', nargs='+',
                    help='The fam file(s) that contains the sample IDs, '
                    'which is an output file from Plink.The files order must '
                    'correspond to the order of Q files')


args = parser.parse_args()

# Define the color palette
if args.color_pallet is None:
    color_pallet = [
        "#648fff",
        "#785ef0",
        "#dc267f",
        "#fe6100",
        "#ffb000",
        "#FFA07A",
        "#FF69B4",
        "#00FF7F",
        "#1E90c8",  # changed from #1E90ff
        "#DC143C",
        "#8B008B",
        "#32CD32",
        "#00CED1",
        "#FF1493",
        "#9ACD32",
        "#FF00FF",
        "#4B0082",
        "#FFD700"]
    color_pallet_version = '# 2024-02-05'
else:
    with open(argparse.color_pallete, 'r') as fi:
        color_pallet = [i.strip() for i in fi.readlines()]
        color_pallet_version = 'Customized'

print(wrap("Using the color pallet version:"+color_pallet_version))


# Read IDs from .fam and id_map files
ids_vcf_full = []
ids_vcf_full_ls = []
fam_count = []
for fam in args.fam:
    with open(fam, 'r') as fi:
        ids_vcf = [line.strip().split()[1] for line in fi.readlines()]
        ids_vcf_full += ids_vcf
        ids_vcf_full_ls.append(ids_vcf)
        fam_count.append(len(ids_vcf))

with open(args.id_map, 'r') as fi:
    id_vcf_uni_map = {}
    id_uni_vcf_map = {}
    lines = fi.read().strip().split('\n')
    if lines[0].startswith('#'):
        lines = lines[1:]
    for line in lines:
        id1, id2 = line.split()
        if id1 in ids_vcf_full:
            id_vcf_uni_map[id1] = id2
            id_uni_vcf_map[id2] = id1
        elif id2 in ids_vcf_full:
            id_vcf_uni_map[id2] = id1
            id_uni_vcf_map[id1] = id2
        else:
            warn(wrap(f'WARNING: Incompatible id_map and fam file'))
            sys.exit(1)

# Read group orders if defined
if args.order_group is not None:
    with open(args.order_group, 'r') as fi:
        order_group = fi.read().strip().split('\n')
else:
    order_group = None

# Read group maps
if args.group_map is not None:
    group_map = {}
    id_group_map = {}
    with open(args.group_map, 'r') as fi:
        lines = fi.read().strip().split('\n')
        if lines[0].startswith('Species') or lines[0].startswith('#'):
            lines = lines[1:]
        for line in lines:
            line = line.split()
            group_map[line[0]] = group_map.get(line[0], [])
            group_map[line[0]].append(line[1])
            id_group_map[line[1]] = line[0]
else:
    group_map = None
    id_group_map = None

# Read Q files
Qs = []
Qk = []
for i, q_file in enumerate(args.Q_file):
    with open(q_file, 'r') as fi:
        Q = fi.read().strip().split('\n')
        k = len(Q[0].split())
        if len(Q) != fam_count[i]:
            raise ValueError('Incompatible Q file and fam file')
        q_dic = {}
        for m in range(fam_count[i]):
            q_dic[ids_vcf_full_ls[i][m]] = Q[m]
        Qs.append(q_dic)
        Qk.append(k)

# Get the output order of samples
samples_order = None
if order_group is not None and group_map is not None:
    samples_order = []
    for group in order_group:
        try:
            ids = group_map[group]
            for id in ids:
                samples_order.append(id)
        except KeyError:
            continue
    for i in Qs:
        if len(i.keys()) > len(samples_order):
            print(len(i.keys(), len(samples_order)))
            raise ValueError('More samples in Q_file than group_order')

if samples_order is None:
    samples_order = set([i for i in ids_vcf_full])

# Output file
header = \
    'ID_vcf\tID_uni\tGroup\t' \
    + "\t".join([f'BestCluster_dataset{str(i+1)}\tColor_dataset{str(i+1)}' 
                 for i in range(len(Qs))]) + "\t"\
    + '\t'.join([f'dataset{str(i+1)}_Cluster{str(k+1)}' 
                 for i in range(len(Qs)) for k in range(Qk[i])])\
    + '\n'\

with open('find_pure_individuals.txt', 'w') as fo:
    fo.write(header)
    for id in samples_order:
        if id in id_vcf_uni_map.keys():
            uni_id = id_vcf_uni_map[id]
            vcf_id = id
        elif id in id_uni_vcf_map.keys():
            uni_id = id
            vcf_id = id_uni_vcf_map[id]
        if id_group_map is not None:
            group = id_group_map[uni_id]
        else:
            raise ValueError('Incompatible group_map and fam file')

        cluster_and_color = ''
        Q_values = ''
        for i, Qi in enumerate(Qs):
            if vcf_id in Qi.keys():
                Q = Qi[vcf_id]
                Q = Q.split()
                Q = [float(i) for i in Q]
                best_cluster = 'Admixed'
                admixed_cluster = ''
                color = '"#888888"'
                for j in range(Qk[i]):
                    if Q[j] >= threshold:
                        best_cluster = 'Cluster'+str(j+1)
                        color = f'"{color_pallet[j]}"'
                        break
                    elif Q[j] > 0:
                        admixed_cluster += '+Cluster'+str(j+1)
                if best_cluster == 'Admixed':
                    best_cluster += admixed_cluster
                if Q_values == '':
                    Q_values = '\t'.join([str(i) for i in Q])
                else:
                    Q_values = Q_values + '\t' + '\t'.join([str(i) for i in Q])

            elif vcf_id not in Qi.keys():
                best_cluster = 'NA'
                color = 'NA'
                if Q_values == '':
                    Q_values = '\t'+'\t'.join(['NA' for i in range(Qk[i])])
                else:
                    Q_values =\
                        Q_values+'\t'+'\t'.join(['NA' for i in range(Qk[i])])
            cluster_and_color += f'{best_cluster}\t{color}\t'


        line = f'{vcf_id}\t{uni_id}\t{group}\t{cluster_and_color}{Q_values}\n'
        fo.write(line)


# # Check ids in fam and id_map
# ## TODO Using enumerate to check the index of ids
# if len(ids) != len(uni_ids) or len(ids) != len(map_ids):
#     print(len(ids), len(uni_ids), len(map_ids))
#     warn(wrap(f'ERROR: Incompatible two input files'))
#     sys.exit(1)
# for id, map_id in zip(ids,map_ids):
#     if id != map_id:
#         warn(wrap(f'WARNING: Incompatible two input files'))

# # Create a dictionary mapping uni_ids to ids
# uni_ids = {id:uni_id for uni_id, id in zip(uni_ids,map_ids)}


# # Find the non-admixture individuals
# fo = open(output_file,'w') 
# fo_Q = open(output_Q,'w')

# with open(input_file_Q, 'r') as fi:
#     n = 0
#     for line in fi:
#         line = line.strip()
#         if line == '':
#             continue

#         non_admix = False
#         n+=1
#         line = line.strip()
#         line = line.replace(' ','\t')
#         Qs = line.split()

#         if n == 1:
#             clusters = "Cluster"+'\tCluster'.join([str(c+1) for c in range(len(Qs))])
#             header=f"#ID\tUni_ID\t{clusters}\n"
#             fo_Q.write(header)
#             header2=f"#ID\tSpecies\tUni_ID\tBestCluster\t{clusters}\n"
#             fo.write(header2)

#         id = ids[n-1]
#         uni_id = uni_ids[id]
#         fo_Q.write(f"{id}\t{uni_id}\t{line}\n")
#         for i in range(len(Qs)):
#             if float(Qs[i]) >= threshold:
#                 non_admix = True
#                 fo.write(f"{id}\t{uni_id[:4]}\t{uni_id}\tCluster{str(i+1)}\t{line}\n")
#                 break
# fo.close()
# fo_Q.close()



# if n != len(ids):
#     print(n, len(ids))
#     warn(wrap('WARNING: Incompatible two input files'))
#     exit


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))