#!/bin/bash

WORKDIR="/work/zruilin/pear/run_FastStructure"

BED="/work/zruilin/pear/run_Population_genetic_filter/output/batch4.remove_b3_remove_Pashia_Pseudopshia/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k"

OUTDIR="$WORKDIR/output/batch4.remove_b3_Pashia_Pseudopshia"

mkdir -p $OUTDIR

perl ./mk_faststructure_genotoul_cv.pl -b $BED -minK 2 -maxK 15 -r 20 -o $OUTDIR

