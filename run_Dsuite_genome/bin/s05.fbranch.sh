#!/bin/bash


#SBATCH -J s05.fbranch.sh
#SBATCH -o s05.fbranch.sh.%J.out
#SBATCH -e s05.fbranch.sh.%J.err


# By Xilong CHEN
# Create date: 2022-10-17
# Contact: chen_xilong@outlook.com

# Adapted by Yuqi.

BIN=$(pwd)
WORK=$(dirname $BIN)
INPUT="$WORK/input"


source $BIN/s00.load_environment.sh

RUNNAME="root_P.cauc_Eastern_tree_Eastern_list.geno20"
OUTPUT2="$WORK/output/$RUNNAME"
OUTPUT5="$WORK/output/s05.fbranch/$RUNNAME"
mkdir -p $OUTPUT5

Dsuite Fbranch \
    $OUTPUT2/root_cauc_Eastern_tree.nwk \
    $OUTPUT2/DTparallel_s02.id_species.Eastern.outgroup_cauc_root_P.cauc_Eastern_tree_Eastern_list.geno20_combined_tree.txt \
    > "${OUTPUT5}/${RUNNAME}_Fbranch.txt"

dtools.py \
    --tree-label-size 10 \
    -n ${OUTPUT5}/${RUNNAME} \
    "${OUTPUT5}/${RUNNAME}_Fbranch.txt" \
    $OUTPUT2/root_cauc_Eastern_tree.nwk


####


RUNNAME="root_P.ussu_Western_tree_Western_list.geno20"
OUTPUT2="$WORK/output/$RUNNAME"
OUTPUT5="$WORK/output/s05.fbranch/$RUNNAME"
mkdir -p $OUTPUT5

TREE=$OUTPUT2/root_ussu_Western.nwk

Dsuite Fbranch \
    $TREE \
    $OUTPUT2/DTparallel_s02.id_species.Western.outgroup_cauc_root_P.ussu_Western_tree_Western_list.geno20_combined_tree.txt \
    > "${OUTPUT5}/${RUNNAME}_Fbranch.txt"

dtools.py \
    --tree-label-size 5 \
    -n ${OUTPUT5}/${RUNNAME} \
    "${OUTPUT5}/${RUNNAME}_Fbranch.txt" \
    $TREE