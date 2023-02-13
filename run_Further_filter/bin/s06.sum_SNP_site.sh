#!/bin/bash


input="/work/zruilin/pear/run_Further_filter/output/s05.filter_concat"
cat $input/*.filtered_variant.count.txt | awk '{sum+=$1}END{print sum}' > sum.snps.count.txt
cat $input/*.filtered_concat.count.txt | awk '{sum+=$1}END{print sum}' > sum.allsites.count.txt


snps=$( cat sum.snps.count.txt )
echo "SNP: $snps"
allsites=$( cat sum.allsites.count.txt )
echo "all sites: $allsites"
# Change invariant_sites=$allsites - $snps
# to this:
invariant_sites=$[$allsites - $snps]
echo "invariant sites: $invariant_sites"
