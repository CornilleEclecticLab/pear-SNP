#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.format_ID.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/10 10:54:46
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
os.makedirs(s03_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Format the ID of pep fasta file
for pep_fa in [pyri_pep_fa, comm_pep_fa]:
    with open(pep_fa, "r") as f:
        basename = os.path.basename(pep_fa)
        prefix = os.path.splitext(basename)[0]
        with open(os.path.join(sub_output_dir, f"{prefix}.formatted.fa"), 'w', encoding='utf-8') as fo:
            for line in f:
                if line.startswith(">") and "Gene=" in line:
                    line = ">"+line.strip().split('Gene=')[1] + "\n"
                fo.write(line)



# Format the ID of gff file
os.makedirs(s04_dir, exist_ok=True)

gene_records = {}
with open(os.path.join(s04_dir, "comm_pyri.gff"), 'w', encoding='utf-8') as fo_concat1, \
        open(os.path.join(s04_dir,"pyri_comm.gff"), 'w', encoding='utf-8') as fo_concat2:
    for gff in [pyri_gff, comm_gff]:
        print(wrap79(f"Processing {gff}"))
        with open(gff, "r") as f:
            basename = os.path.basename(gff)
            prefix = os.path.splitext(basename)[0]
            pre_prefix = prefix.split('.')[0]
        
            chr_map = {}
            chr_map_file = os.path.join(input_dir, f"{pre_prefix}.chrID_map.txt")
            if os.path.exists(chr_map_file):
                with open(chr_map_file, 'r') as map:
                    for line in map:
                        lines = line.strip().split('\t')
                        chr_map[lines[0]] = lines[1]
        
            for line in f:
                if line.startswith("#") or line.strip() == "":
                    continue 
                
                lines = line.strip().split('\t')
                chr = lines[0]
                start = int(lines[3]) - 1
                end = int(lines[4])
                annotation = lines[8]
                if chr.startswith("Chr"):
                    chr = chr.replace("Chr", "co")
                    gene = annotation.split('Name=')[1].split(';')[0]
                elif chr in chr_map:
                    chr = chr_map[chr].replace("Chr", "py")
                    gene = annotation.split('Accession=')[1].split(';')[0]
                else:
                    raise ValueError(f"Unknown chromosome ID: {chr}")

                if gene in gene_records:
                    print(f"Warning: Skiped duplicate gene ID: {gene}")
                    continue
                else:
                    gene_records[gene] = True
                    fo_concat1.write(f"{chr}\t{gene}\t{start}\t{end}\n")
                    fo_concat2.write(f"{chr}\t{gene}\t{start}\t{end}\n")
        print(wrap79(f"Processed {prefix}.bed with {len(gene_records)} accumulated unique genes"))
fo_concat1.close()
fo_concat2.close()

end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
