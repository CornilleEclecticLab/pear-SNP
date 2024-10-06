#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : phy_invariant_cleaner.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/10/03 11:46:14
# @Description:
#    This script removes invariant sites from a PHYLIP file and writes only variants to a new file,
#    to suppress the rejection of IQ-TREE (+ASC models) due to invariant sites.

#    - Compatible with vcf2phylip.py v2.9 and PHYLIP format (Relaxed and Nucleic acid only).
#    - vcf2phylip.py: https://github.com/edgardomortiz/vcf2phylip/
#    - Phylip format: https://www.phylo.org/index.php/help/phylip
#    - IUPAC code: https://www.phylo.org/index.php/help/iupac


import argparse
import textwrap
import datetime
from Bio import AlignIO

print(f'{" Start ":=^79}')

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' ' * 4)


# Get arguments
parser = argparse.ArgumentParser(description="Remove invariant sites from a PHYLIP file.")
parser.add_argument("-i", "--input", type=str, required=True, help="Input PHYLIP filename")
parser.add_argument("-o", "--outfile", type=str, default="phy_invariant_cleaner.variants.phy", help="Output filename")
parser.add_argument("-n", "--no-check", action="store_true", help="Do not check for unknown bases")
args = parser.parse_args()

# Start time
start_time = datetime.datetime.now()

# Read the input PHYLIP alignment
input_file = args.input
with open(input_file) as fi:
    alignment = AlignIO.read(fi, "phylip-relaxed")
    alignment_array = [list(rec) for rec in alignment]
    alignment_transposed = list(zip(*alignment_array))


# Define the general bases and IUPAC ambiguous bases
GENERAL_BASES = ["A", "C", "G", "T", "U", 
                 "R", "Y","K", "M", "S", "W", "B", "D", "H", "V", "N"]
IUPAC_AMBIGUOUS_BASES = ["R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N"]
IUPAC_AMBIGUOUS_DICT= {
    'R':'GA',  'Y':'TC',  'K':'GT',  'M':'AC', 'S':'GC', 'W':'AT',
    'B':'GTC', 'D':'GAT', 'H':'ACT', 'V':'GCA','N':'AGCT'} 


# Function to check for unknown bases
def check_unknown_bases(bases_set, position):
    if args.no_check:
        return
    for base in bases_set:
        if base not in GENERAL_BASES and base not in ["-", "N", "X"]:
            raise ValueError(
                wrap79(f"Unknown base '{base}' in the alignment at position {position + 1}. The accepted bases are {set(GENERAL_BASES+['-','N','X'])} .")
            )


# Filter out invariant sites for polymorphic sites
# invariant_sites_column = []
variant_sites = []

for i, codes in enumerate(alignment_transposed):
    # Set of unique bases in the site, ignoring gaps/missing data ('-', 'N', 'X')
    unique_bases = set([c.upper() for c in codes]).difference(set(['-', 'N','X']))
    
    if len(unique_bases) == 0:
        # Treat as invariant since it contains no true information
        # invariant_sites_column.append(i)
        pass
    
    elif len(unique_bases) == 1:           
        # Check for unknown bases
        check_unknown_bases(unique_bases, i)
        # If there's only one unique base, the site is invariant
        # invariant_sites_column.append(i)
        pass
    
    elif len(unique_bases) == 2:
        # Check for unknown bases
        check_unknown_bases(unique_bases, i)
        # Handle IUPAC ambiguous bases
        hits_ambiguous = [base for base in unique_bases if base in IUPAC_AMBIGUOUS_BASES]
        hits_non_ambiguous = [base for base in unique_bases if base not in IUPAC_AMBIGUOUS_BASES]  # non-ambiogous = 'A' 'T' 'G' 'C' 'U'
        if len(hits_ambiguous) == 1 and len(hits_non_ambiguous) == 1:
            if hits_ambiguous[0] in IUPAC_AMBIGUOUS_BASES[hits_non_ambiguous[0]]:
                # invariant_sites_column.append(i)
                pass
            else:
                variant_sites.append(codes)
        else:
            variant_sites.append(codes)
    else:
        variant_sites.append(codes)


# Transpose back to the original orientation
filtered_alignment = list(zip(*variant_sites))


# Write to a new phylip file with variant-only sites
output = args.outfile
with open(output, "w") as f:
    f.write(f"{len(alignment)} {len(variant_sites)}\n")
    for i, record in enumerate(alignment):
        seq_id = record.id
        sequence = ''.join(filtered_alignment[i])
        f.write("{id:<10} {seq}\n".format(id=seq_id, seq=sequence))


# End time
end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str('Elapsed time: '+str(end_time-start_time)).center(79))
