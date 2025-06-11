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

gff_file_dic = {
    'PPY': os.path.join(input_dir, 'GWHBAOS00000000.no_blank.chr.gff'),
    'PCOM': os.path.join(input_dir, 'PyrusCommunis_BartlettDHv2.0.chr.gff')
}

mask_pass_bed_dic ={
    'PPY': os.path.join(input_dir, 'GWHBAOS00000000.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed'),
    'PCOM': os.path.join(input_dir, '../input/PyrusCommunis_BartlettDHv2.0.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed')
}

cds_masked_bed_dic = {
    'PPY': os.path.join(input_dir, 'PPY.cds_masked.bed'),
    'PCOM': os.path.join(input_dir, 'PCOM.cds_masked.bed')
}
