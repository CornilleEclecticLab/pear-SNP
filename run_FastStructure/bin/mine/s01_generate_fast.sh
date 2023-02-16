#!/bin/bash

WORKDIR="/work/zruilin/RNAseq_fast_structure"

BED="$WORKDIR/input/Combine_Chr.geno20_maf005.anno.syno.thin5k"
OUTDIR="$WORKDIR/output/s01_generate_fast"
mkdir $OUTDIR

perl ./mk_fastStructure_genotoul_cv.pl -b $BED -minK 2 -maxK 15 -r 20 -o $OUTDIR

