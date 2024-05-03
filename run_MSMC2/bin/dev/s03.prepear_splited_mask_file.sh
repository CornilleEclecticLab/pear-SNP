#!/usr/bin/env bash

mkdir ../input/mask && cd ../input/mask && \
for i in $(cat ../s02.chromosome_list.txt) ; do echo $i | xargs -I % grep % GWHBAOS00000000.genome.genmap.mask.bed.pass.bed > GWHBAOS00000000.genome.$i.chromosome.genmap.mask.bed.pass.bed ; done

