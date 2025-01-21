#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.count_genes_by_windows_from_gff.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/01/14 18:08:46
# @Description:
#    

import pandas as pd
from intervaltree import IntervalTree
import datetime
import sys
import textwrap
import os
import glob
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "0.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(output_dir,script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
# os.makedirs(sub_output_dir, exist_ok=True)
# os.makedirs(sub_script_dir, exist_ok=True)


# Define window size and step size for counting genes
# Ref https://doi.org/10.1093/g3journal/jkae003
window_size = 1_000_000  # 1 Mb
step_size =  100_000     # 100 kb


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)

# Function to read gff file
def read_gff(gff_file):
    # Load GFF file
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
    df.columns = ["seqname", "source", "feature", "start",
                  "end", "score", "strand", "frame", "attribute"]

    # Filter for genes
    genes = df[df['feature'] == 'gene']
    return genes

# Function to count genes in each window
def count(genes, chromosomes, window_size, step_size):
    
    # Calculate gene counts for each window
    results = []
    for seqname in genes['seqname'].unique():
        if seqname not in chromosomes:
            continue
        seq_length = max(genes[genes['seqname'] == seqname]['end'])
        
        # Filter genes on the same chromosome
        seq_genes = genes[genes['seqname'] == seqname]
        
        # Create Interval Tree for gene coordinates
        tree = IntervalTree()
        
        for start, end in zip(seq_genes['start'], seq_genes['end']):
            tree[start:end] = 1
        
        for start in range(1, seq_length+1, step_size):
            end = start + window_size
            interval_count = len(tree[start:end])
            results.append((seqname, start-1, end-1, interval_count))

    # Convert results to DataFrame
    result_df = pd.DataFrame(
        results, columns=["Chr", "Start", "End", "Gene_count"])

    return result_df

# Function to read candidate chromosomes
def read_chromosomes(chrom_file):
    with open(chrom_file) as f:
        chromosomes = [line.strip() for line in f]
    return chromosomes


# Main
gff_files = glob.glob(os.path.join(input_dir, '*.chr_list.gff'))
for gff_file in gff_files:
    print(wrap79(f'Processing {gff_file} ...'))
    
    basename = os.path.basename(gff_file)
    chrom_file = os.path.join(input_dir, basename.replace(
        '.chr_list.gff', '.chr_list.txt'))
    chromosomes = read_chromosomes(chrom_file)
    
    genes_df = read_gff(gff_file)
    
    result_df = count(genes_df, chromosomes, window_size, step_size)
    
    output_file = os.path.join(output_dir, os.path.basename(gff_file).replace(
        '.chr_list.gff', f'.wind{window_size}_step{step_size}.tsv'))
    result_df.to_csv(output_file, sep='\t', index=False)



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
