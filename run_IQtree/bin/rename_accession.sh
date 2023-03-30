#!/bin/bash

cp "$1" "$1".rename

while read key value; do
    sed -i "s/\b$key\b/$value/g" "$1".rename
done < accession_to_id.txt

# while read key value; do
#     sed  "s/\b$key\b/$value/g" $1 > $i.rename
# done < accession_to_id.txt 
