#!/usr/bin/bash

source s00.load_bedtools.sh


bedtools slop \
    -i ../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.gff \
    -g ../input/karyotype.txt \
    -b 10000 \
    -s  \
> ../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop10kb_up_down.gff
