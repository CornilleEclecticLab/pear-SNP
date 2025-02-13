#!/usr/bin/bash

source s00.load_bedtools.sh


bedtools slop \
    -i ../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.gff \
    -g ../input/karyotype.txt \
    -l 1000 \
    -r 0 \
    -s  \
> ../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop1kb.gff 
