#!/bin/bash

WORKDIR="/work/zruilin/pear/run_FastStructure"

BED="/work/zruilin/pear/run_Population_genetic_filter/output/batch2.removeCloneAnd3LowQualEuropeWildPear/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k"
OUTDIR="$WORKDIR/output/batch2.removeCloneAnd3LowQualEuropeWildPear"

mkdir -p $OUTDIR

perl ./mk_faststructure_genotoul_cv.pl -b $BED -minK 2 -maxK 15 -r 20 -o $OUTDIR

