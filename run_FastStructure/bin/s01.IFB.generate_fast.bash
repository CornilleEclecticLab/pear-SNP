#!/bin/bash

WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_FastStructure"


BED_PREFIX="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/branch8.batch15.pear_Jul2023_noclone/pear_Jul2023_noclone.Combine_Chr.geno20_maf005.anno.syno.thin8k"

OUTDIR="$WORKDIR/output/branch8.bathch15.pear_Jul2023_noclone"



perl ./mk_faststructure_IFB_cv.pl -b $BED_PREFIX -minK 2 -maxK 15 -r 30 -o $OUTDIR

