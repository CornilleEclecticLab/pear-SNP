#!/bin/bash
# By Xilong CHEN
# Create date: 2022-08-22
# Contact: chen_xilong@outlook.com
# Adaptated by Yuqi NIE on 2023-11-29 17:37:08

#SBATCH -J s03.combine
#SBATCH -o s03.combine.NOGIT.%J.out
#SBATCH -e s03.combine.NOGIT.%J.err

BIN=$(pwd)
WORK=$(dirname $BIN)
INPUT="$WORK/input"
OUTPUT="$WORK/output/s03.combine"

source $BIN/s00.load_environment.sh

Dsuite DtriosCombine \
-o $OUTPUT/genome \
-t $INPUT/tree.txt \