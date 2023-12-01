#!/bin/bash
# By Xilong CHEN
# Create date: 2022-10-17
# Contact: chen_xilong@outlook.com

#

BIN=$(pwd)
WORK=$(dirname $BIN)
INPUT="$WORK/input"
OUTPUT="$WORK/output/s05.fbranch"
OUTPUT3="$WORK/output/s03.combine"

source $BIN/s00.load_environment.sh

Dsuite Fbranch \
    $INPUT/tree_astral.txt \
    $OUTPUT3/genome_combined_tree.txt \
    > $OUTPUT/tree_astral_with_geneflow_Fbranch.txt

dtools.py \
    -n $OUTPUT/tree_astral_with_geneflow_Fbranch \
    $OUTPUT/tree_astral_with_geneflow_Fbranch.txt \
    $INPUT/tree_astral.txt


Dsuite Fbranch \
    $INPUT/tree_twisst.txt \
    $OUTPUT3/genome_combined_tree.txt \
    > $OUTPUT/tree_twisst_with_geneflow_Fbranch.txt

dtools.py \
    -n $OUTPUT/tree_twisst_with_geneflow_Fbranch \ 
    $OUTPUT/tree_twisst_with_geneflow_Fbranch.txt \
    $INPUT/tree_twisst.txt

