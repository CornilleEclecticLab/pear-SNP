#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s08.count_collinearity_genes_num.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/12 15:45:12
# @Description:
#    

import datetime
import sys
import textwrap
import os
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

# Read collinearity file
with open(os.path.join(s05_dir, "intersection.collinearity.txt"), 'r') as f:
    lines = f.readlines()
    intersection_genes1= {line.strip().split()[0]: True for line in lines}
    intersection_genes2= {line.strip().split()[1]: True for line in lines}
    intersection_genes_dual = [intersection_genes1, intersection_genes2]
    

with open(os.path.join(s05_dir, "union.collinearity.txt"), 'r') as f:
    lines = f.readlines()
    union_genes1 = {line.strip().split()[0]: True for line in lines}
    union_genes2 = {line.strip().split()[1]: True for line in lines}
    union_genes_dual = [union_genes1, union_genes2]


# Read gff
for gff in [comm_gff, pyri_gff]:
    with open(gff, 'r') as f:
        gff_genes = {}
        basename = os.path.basename(gff)
        prefix = '.'.join(str(basename).split('.')[0:-3])
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            elif line == '':
                continue
            elif line.split('\t')[2] != 'gene':
                continue
            annotation = line.split('\t')[8]
            annotation_parse = { key_value.split('=')[0]: key_value.split('=')[1] for key_value in annotation.split(';') if '=' in key_value}
            if 'Name' in annotation_parse:
                gene_name = annotation_parse['Name']
            elif 'Accession' in annotation_parse:
                gene_name = annotation_parse['Accession']
            elif 'ID' in annotation_parse:
                gene_name = annotation_parse['ID']
            else:
                raise ValueError(f"Unknown gene name fromat found in {annotation}")
            gff_genes[gene_name] = True
                
    # Count total genes
    print(wrap79(f"{len(gff_genes)}: total genes in {basename}. "))
    
    for i in range(len(intersection_genes_dual)):
        if all (gene in gff_genes for gene in intersection_genes_dual[i]):
            intersection_genes_num = len(intersection_genes_dual[i])
            print(wrap79(f"{intersection_genes_num}: Intersection collinearity genes in {basename}"))

            with open(os.path.join(sub_output_dir, f"{prefix}.intersection.collinearity.txt"), 'w') as f_c,\
                open(os.path.join(sub_output_dir, f"{prefix}.intersection.non_collinearity.txt"), 'w') as f_n,\
                open(os.path.join(sub_output_dir, f"{prefix}.gff.txt"), 'w') as f_g:
                    for gene in gff_genes:
                        f_g.write(f"{gene}\n")
                        if gene in intersection_genes_dual[i]:
                            f_c.write(f"{gene}\n")
                        else:
                            f_n.write(f"{gene}\n")
        
    
    for i in range(len(union_genes_dual)):
        if all (gene in gff_genes for gene in union_genes_dual[i]):
            union_genes_num = len(union_genes_dual[i])
            print(wrap79(f"{union_genes_num}: Union collinearity genes in {basename}"))

            with open(os.path.join(sub_output_dir, f"{prefix}.union.collinearity.txt"), 'w') as f_c,\
                    open(os.path.join(sub_output_dir, f"{prefix}.union.non_collinearity.txt"), 'w') as f_n:
                    for gene in gff_genes:
                        if gene in union_genes_dual[i]:
                            f_c.write(f"{gene}\n")
                        else:
                            f_n.write(f"{gene}\n")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))