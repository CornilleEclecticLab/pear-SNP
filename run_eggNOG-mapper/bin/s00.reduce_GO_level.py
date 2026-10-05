#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s00.reduce_GO_level.py
# @Time(CET): 2025/04/29 14:25:46
# @Description:
#    This script is modified from Somia's
#    To load "go-basic.obo" file : https://current.geneontology.org/ontology/go-basic.obo

import goatools
import pandas as pd
from goatools.obo_parser import GODag


import datetime
import sys
import textwrap
import os
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


def get_leveln_terms(go_list_str, go_dag, nLevel):

	if pd.isna(go_list_str) or go_list_str == "None":
		return "None"

	go_ids = [go.strip() for go in go_list_str.split(',')]
	leveln_terms = set()

	for go_id in go_ids:
		try:
			term = go_dag[go_id]
			ancestors = term.get_all_parents()
			for ancestor_id in ancestors:
				ancestor_term = go_dag[ancestor_id]
				if ancestor_term.depth == nLevel:
					leveln_terms.add(ancestor_term.name)
		except KeyError:
			continue

	return ",".join(sorted(leveln_terms)) if leveln_terms else "None"


if __name__ == "__main__":

	# Load GO database and GO Directed Acyclic Graph
	print("\ngoatools version : ", goatools.__version__)
	go_obo_file = "../input/go-basic.obo"
	go_dag = GODag(go_obo_file)

	# Example
	data = {
		'GOs': [
			'GO:0003674,GO:0005515,GO:0007049',
			'GO:0003735,GO:0006412',
			'GO:0005634',
			'-'
                ]
	}
	df = pd.DataFrame(data)
	print("\n", df)

	# Convert
	df['GOs_level0'] = df['GOs'].apply(
	    lambda gos: get_leveln_terms(gos, go_dag, 0))
	df['GOs_level1'] = df['GOs'].apply(
	    lambda gos: get_leveln_terms(gos, go_dag, 1))
	df['GOs_level2'] = df['GOs'].apply(
	    lambda gos: get_leveln_terms(gos, go_dag, 2))
	df['GOs_level3'] = df['GOs'].apply(
	    lambda gos: get_leveln_terms(gos, go_dag, 3))

	# df.to_csv(f + ".GOs", sep="\t", index=False)
	print("\n", df)
 
	# Example 2
	data2 = 'GO:0003674,GO:0005515,GO:0007049'
	print(get_leveln_terms(data2, go_dag, 0))
	print(get_leveln_terms(data2, go_dag, 1))
	print(get_leveln_terms(data2, go_dag, 2))
	print(get_leveln_terms(data2, go_dag, 3))


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))