#!/usr/bin/env bash

# By Yuqi NIE, 2023-07-17

# This script is used to record the job ID of each job submitted to the cluster.
for i in *.out; do
	base=$(echo "$i" | awk -F '.' '{for (n=1; n<=(NF-2); n++) printf "%s.", $n; printf "sh"}');
	id=$(echo "$i" | awk -F '.' '{print $(NF-1)}');
	printf "%s\t%s\n" "$base" "$id";
done > job_ID.txt