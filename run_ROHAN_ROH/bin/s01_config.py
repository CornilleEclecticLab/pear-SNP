#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01_config.py
# @Author   : NIE Yuqi
# @Time(Europe/Paris): 2026/04/30 10:31:00
# @Description:
#    

import datetime
import sys
import textwrap
import os

version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(output_dir,script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)


sample_population_file = os.path.join(input_dir, 's01.sample_tab_population.txt')
bed_file = os.path.join(input_dir, 'GWHBAOS00000000.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed')
vcf_file_list = os.path.join(input_dir, 'vcf_file_list.txt')
cram_file_list = os.path.join(input_dir, 'cram_file_list.txt')
ref_genome = os.path.join(input_dir, 'GWHBAOS00000000.genome.fasta')
chr_list_file = os.path.join(input_dir, 's01.chr_list.txt')
tstv_dir = os.path.join(output_dir, 's01.calculate_tstv_from_vcf')



# Read sample population file
sample_population_dic = dict()
population_sample_dic = dict()
with open(sample_population_file, 'r') as fi:
    for line in fi:
        line = line.strip().split()
        s, p = line[0], line[1]
        if s not in sample_population_dic:
            sample_population_dic[s] = p
        else:
            raise ValueError(f'Dup sample {s}')
        
        if p not in population_sample_dic:
            population_sample_dic[p] = []
        population_sample_dic[p].append(s)

# Read chrom list 
chr_list = []
with open(chr_list_file, 'r') as fi:
    chr_list = [line.strip() for line in fi]

