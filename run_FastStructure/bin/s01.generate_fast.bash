#!/bin/bash

WORKDIR="/work/zruilin/Peach/run_FastStructure"

BED="$WORKDIR/input/Combine_Chr.geno20_maf005.anno.syno.thin8k"
OUTDIR="$WORKDIR/output/s01.generate_fast"
mkdir $OUTDIR

perl ./mk_faststructure_genotoul_cv.pl -b $BED -minK 2 -maxK 10 -r 20 -o $OUTDIR

