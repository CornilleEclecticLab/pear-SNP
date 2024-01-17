#!/bin/bash

source s00.config.sh

mkdir ${OUTPUT_DIR}

perl ./mk_faststructure_IFB_cv.pl \
    -b $BED_PREFIX \
    -minK $MIN_K \
    -maxK $MAX_K \
    -r 30 \
    -o ${OUTPUT_DIR}
