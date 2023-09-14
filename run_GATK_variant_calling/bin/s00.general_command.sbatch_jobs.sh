#!/usr/bin/env bash

# By Yuqi NIE, on 2023-07-17

# for i in *.sh ; do sbatch -A pear_snp2 -p long -c 10 --mem=28G $i | awk '{print "'$i'\t"$4}' ; done > job_ID.txt
