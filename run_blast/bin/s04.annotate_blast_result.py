#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.annotate_blast_result.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/15 17:51:03
# @Description:
#    

import datetime
import sys
import textwrap
import os
import glob
from s00_config import *

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(output_dir,script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)
# os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)



# Loading db annotation
db_annotation = {} # {sp: {gene_id: [symbol, description]}}
for db_file in glob.glob(f'{s03_dir}/*.annotation.tsv'):
    with open(db_file, 'r') as fi:
        sp = os.path.basename(db_file).split('.')[0]
        db_annotation[sp] = {}
        for line in fi:
            gene_id, symbol, description = line.strip().split('\t')
            db_annotation[sp][gene_id] = [symbol, description]


def parse_blast(line):
    lines = line.strip().split()
    query_id = lines[0]
    subject_id = lines[1]
    return query_id, subject_id

blast_files = glob.glob(f'{s02_dir}/*.blastp.tsv')
blast_dic = {} # {query_sp: {qeury_gene: {target_sp: target_gene}}}

                  
for blast_file in blast_files:
    basename = os.path.basename(blast_file)
    query_sp = '.'.join(basename.split('.')[:-3])
    target_sp = basename.split('.')[-3]

    print(f'Reading {basename} {query_sp} {target_sp} ...')
    
    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip()
            query_id, subject_id = parse_blast(line)
            blast_dic[query_sp] = blast_dic.get(query_sp, {})
            blast_dic[query_sp][query_id] = blast_dic[query_sp].get(query_id, {})
            blast_dic[query_sp][query_id][target_sp] = subject_id
    

for query_sp, query_sp_values in blast_dic.items():
    target_sps = set(target_sp
                     for query_sp_values_value in query_sp_values.values()
                     for target_sp in query_sp_values_value.keys())
    target_sps = sorted(target_sps)
    
    with open(os.path.join(sub_output_dir, f'{query_sp}.blastp.annotated.tsv'), 'w') as fo:
        header = 'Query\t' + '\t'.join([f'{sp}_best_target\t{sp}_symbol\t{sp}_description' for sp in target_sps]) + '\n'
        fo.write(header)

        for query_gene, target_sp_dic in query_sp_values.items():
            fo.write(f'{query_gene}\t')
            for target_sp in target_sps:
                target_gene = target_sp_dic.get(target_sp, None)
                if target_gene is None:
                    fo.write("-\t-\t-\t")
                else:
                    [symbol, description] = db_annotation[target_sp][target_gene]
                    fo.write(f'{target_gene}\t{symbol}\t{description}\t')
            fo.write('\n')

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))