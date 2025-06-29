#!/bin/bash

bcftools view <vcf> \
| bcftools +fill-tags -- -t F_MISSING \
| bcftools query -f '%CHROM\t%POS\t%INFO/F_MISSING\n' > fmiss.txt
