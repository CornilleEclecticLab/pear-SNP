#!/usr/bin/env python3
import os
import sys

script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)


# This table header is sample_ID_in_VCF, sampleID_Uniform, population, depth, reads_length
individuals_table = os.path.join(input_dir, "s01.individuals_table.txt") 

# This directory contains the CRAM files, check the reference genome and the index files
# cram file name expected format: *.sampleID.*.*.cram
cram_dir = os.path.join(input_dir, "cram_files", "*.cram")
