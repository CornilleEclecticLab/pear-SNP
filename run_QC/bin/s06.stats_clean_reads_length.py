#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.stats_clean_reads_length.py
# @Version  : 1.1.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/04/10 15:12:15
# @Description:
#     Update: 1.1.0:
#             support multi QC v1.18

import datetime
import sys
import textwrap
import os
import pandas as pd
import glob



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



def read_table(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df

file_paths = glob.glob("../input/s06.input/s06_input_data/multiqc_general_stats.txt")


output_file = os.path.join(output_dir, "s06.stats_clean_reads_length.txt")
output_file2 = os.path.join(output_dir, "s06.stats_clean_reads_pair_average_length.txt")

with open(output_file, "w") as f:
    with open(output_file2, "w") as f2:
        for file_path in file_paths:
            table = read_table(file_path)

            # 1 Remove leading spaces and "clean." from Sample IDs
            table['Sample'] = table['Sample'].str.lstrip().str.replace('^clean\.', '', regex=True).str.replace("FKDN23.+?_L1_","",regex=True)

            # 2 Compare lengths for pairs of Sample IDs ending with "_1" and "_2"
            table['_1'] = table['Sample'].str.endswith('_1')
            table['_2'] = table['Sample'].str.endswith('_2')
            sample_pairs = table[table['_1'] | table['_2']].groupby(table['Sample'].str.replace(r'(_1|_2)$', '', regex=True))
            for name, stats in sample_pairs:
                if len(stats) == 2:
                    length_1 = stats.loc[stats['_1'], 'avg_sequence_length'].values[0]
                    length_2 = stats.loc[stats['_2'], 'avg_sequence_length'].values[0]

                    if abs(length_1 - length_2) > 1:
                        print(f"Sample pair {name}_1 and {name}_2 has length difference > 1\n")

                    avg_length = (length_1 + length_2) / 2
                    f2.write(f"{name}\t{avg_length}\n")

            # 3 Check for duplicate Sample IDs with inconsistent lengths
            duplicate_samples = table[table['Sample'].duplicated(keep=False)]
            for name, stats in duplicate_samples.groupby('Sample'):
                if len(stats['avg_sequence_length'].unique()) > 1:
                    print(f"Sample {name} has inconsistent lengths: {', '.join(map(str, stats['avg_sequence_length'].unique()))}\n")


            selected_columns = table[["Sample","avg_sequence_length"]]
            f.write(selected_columns.to_string(index=False) + "\n\n")
