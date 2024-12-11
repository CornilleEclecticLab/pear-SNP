#!/bin/bash

OUT=../input/s01.6pops.3sample.outgroup_cauc.txt
SAMPLES=../input/sample_tab_population.txt
POPS=../input/pops

if [ -f $OUT ]; then
    echo "sampling results file already exists."
    exit
fi

while read -r pop; do
  grep  "${pop}$" $SAMPLES | shuf -n 3 >> $OUT
done < $POPS
