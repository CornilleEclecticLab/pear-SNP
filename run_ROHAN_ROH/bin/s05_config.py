#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s05_config.py
# @Author   : NIE Yuqi
# @Description:
#   Config for s05_mutation_burden_by_ROH.py.

import os
import sys

config_path  = os.path.dirname(os.path.realpath(__file__))   # .../run_ROHAN_ROH/bin
ROHAN_DIR    = os.path.dirname(config_path)                  # .../run_ROHAN_ROH
ROHAN_input  = os.path.join(ROHAN_DIR, 'input')
ROHAN_output = os.path.join(ROHAN_DIR, 'output')

ROHAN_DIR_PCOM    = os.path.join(os.path.dirname(ROHAN_DIR), 'run_ROHAN_ROH_refer_on_communis/')                  # .../run_ROHAN_ROH
ROHAN_input_PCOM  = os.path.join(ROHAN_DIR_PCOM, 'input')
ROHAN_output_PCOM = os.path.join(ROHAN_DIR_PCOM, 'output')

SIFT4G_DIR    = os.path.join(os.path.dirname(ROHAN_DIR), 'run_SIFT4G')
SIFT4G_input  = os.path.join(SIFT4G_DIR, 'input')
SIFT4G_output = os.path.join(SIFT4G_DIR, 'output')

ref_genomes = ['PPY', 'PCOM']

# inputs from run_SIFT4G
chr_list_file_dic = {
    'PPY':  os.path.join(SIFT4G_input, 'GWHBAOS00000000.chr_list.txt'),
    'PCOM': os.path.join(SIFT4G_input, 'PyrusCommunis_BartlettDHv2.0.chr_list.txt'),
}

vcf_file_dir_dic = {
    'PPY':  os.path.join(SIFT4G_output, 's03.get_vcf_by_pop_from_variant_PPY'),
    'PCOM': os.path.join(SIFT4G_output, 's03.get_vcf_by_pop_from_variant_PCOM'),
}

sift4g_dir_dic = {
    'PPY':  os.path.join(SIFT4G_output, 's02_annot', 'PPY'),
    'PCOM': os.path.join(SIFT4G_output, 's02_annot', 'PCOM'),
}

individuals_map_file_dic = {
    'PPY':  os.path.join(SIFT4G_input, 's01.sample_tab_population_PPY.txt'),
    'PCOM': os.path.join(SIFT4G_input, 's01.sample_tab_population_PCOM.txt'),
}

# inputs from run_ROHAN_ROH
roh_summary_dir_dic = {
    'PPY':  os.path.join(ROHAN_output, 's02.generate_run_rohan'),
    'PCOM': os.path.join(ROHAN_output_PCOM, 's02.generate_run_rohan'),
}

def require_paths():
    for ref in ref_genomes:
        for label, d in (('sift4g_dir', sift4g_dir_dic[ref]),
                         ('vcf_file_dir', vcf_file_dir_dic[ref]),
                         ('roh_summary_dir', roh_summary_dir_dic[ref])):
            if not os.path.isdir(d):
                sys.exit(f'[s05_config] {label} for {ref} not found: {d}')