#!/usr/bin/env python3

import os
import sys

script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')
sub_output_dir = os.path.join(output_dir, script_basename)
sub_script_dir = os.path.join(bin_dir, script_basename)

ref_genomes = ['PPY', 'PCOM']


chr_list_file_dic = {
    'PPY': os.path.join(input_dir, 'GWHBAOS00000000.chr_list.txt'),
    'PCOM': os.path.join(input_dir, 'PyrusCommunis_BartlettDHv2.0.chr_list.txt')
}

vcf_file_dir_dic = {
    'PPY': os.path.join(output_dir, 's03.get_vcf_by_pop_from_variant_PPY'),
    'PCOM': os.path.join(output_dir, 's03.get_vcf_by_pop_from_variant_PCOM')
}

sift4g_dir_dic = {
    'PPY': os.path.join(output_dir, 's02_annot', 'PPY'),
    'PCOM': os.path.join(output_dir, 's02_annot', 'PCOM')
}

individuals_map_file_dic = {
    'PPY': os.path.join(input_dir, 's01.sample_tab_population_PPY.txt'),
    'PCOM': os.path.join(input_dir, 's01.sample_tab_population_PCOM.txt')
}
