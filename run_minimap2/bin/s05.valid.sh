#!/usr/bin/bash

echo "Minimap2 results md5:"
head -n 17 ../output/aln.id.count.txt | awk '{print $2"\t"$1}' | sort | sed "s/Chr//g" | md5sum

echo "Fasta annotations md5:"
cat ../input/chrID_number.txt | sort | md5sum
