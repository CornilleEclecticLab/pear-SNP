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


variant_vcf_list_file_dic = {
    'PPY': os.path.join(input_dir, 's03.variant_vcf_list.PPY.txt'),
    'PCOM': os.path.join(input_dir, 's03.variant_vcf_list.PCOM.txt')
}
