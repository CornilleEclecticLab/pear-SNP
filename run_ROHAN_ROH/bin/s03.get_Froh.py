#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.get_roh.py
# @Author   : NIE Yuqi
# @Time(Europe/Paris): 2026/06/05 16:25:23
# @Description:
#    

import datetime
import sys
import textwrap
import os
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.0.0"

import s01_config as cfg
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



# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


import sys
import re
from pathlib import Path

def extract_froh_matrix(roh_base_dir, output_file):
    """
    Glob ROHan summary files, extract F_ROH from each file and output a matrix.
    """
    # ROHan summary path : ../output/s02.*/<pop_name>/<sample_ID>.summary.txt
    search_pattern = "s02.generate_run_rohan/*/*.summary.txt"
    target_files = list(Path(roh_base_dir).glob(search_pattern))
    
    print(f"Found {len(target_files)} summary files.")

    # eg. "Segments in ROH(%) : 0.16129 (0,0.161551)"

    # Genome-wide theta outside ROH: 0.00476966 0.00475731 0.00498224

    # Genome-wide theta inc. ROH: 0.00476623 0.00476389 0.00498169

    # Segments unclassified : 0 (0,0)

    # Segments unclassified (%): 0 (0,0)

    # Segments in ROH : 1000000 (0,1000000)

    # Segments in ROH(%) : 0.16129 (0,0.161551)

    # Segments in non-ROH : 619000000 (619000000,620000000)

    # Segments in non-ROH (%) : 99.8387 (99.8387,100)

    # Avg. length of ROH : 1e+06 (0,1e+06)

    roh_pattern = re.compile(
        r"^Segments in ROH\(%\)\s*:\s*([0-9.eE+-]+)\s*\(\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\)"
    )

    #
    with open(output_file, 'w') as f_out:
        # Header
        f_out.write("Sample\tPopulation\tF_ROH\tF_ROH_lower\tF_ROH_upper\n")
        
        for filepath in target_files:
            sample_name = filepath.name.replace(".summary.txt", "")
            pop_name = filepath.parent.name
            
            froh_coef = "NA"
            froh_lower = "NA"
            froh_upper = "NA"
            
            try:
                with open(filepath, 'r') as f_in:
                    for line in f_in:
                        match = roh_pattern.search(line)
                        if match:
                            froh_percent = float(match.group(1))
                            # Percent to coefficient
                            froh_coef = froh_percent / 100.0
                            
                            lower_percent = float(match.group(2))
                            froh_lower = lower_percent / 100.0
                            
                            upper_percent = float(match.group(3))
                            froh_upper = upper_percent / 100.0
                            break
            except Exception as e:
                print(f"Warning: Failed to parse {filepath}. Reason: {e}")
                
            f_out.write(f"{sample_name}\t{pop_name}\t{froh_coef}\t{froh_lower}\t{froh_upper}\n")

    print(f"\n>>> Matrix saved to {output_file}")

if __name__ == "__main__":
    BASE_DIRECTORY = "../output" 
    OUTPUT_MATRIX = "../output/FrohMatrix.tsv"
    
    extract_froh_matrix(BASE_DIRECTORY, OUTPUT_MATRIX)

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))