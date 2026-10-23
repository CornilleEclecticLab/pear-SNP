#!/bin/bash

# @File     :   s04.Jubail.bootstrap.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi modified from Xilong's scripts
# @Time(CET):   2025/11/20 09:36:34
#Description:
    # 

# Load modules easySFS
# . s00.load_easySFS.sh

set -euo pipefail

### Configure
SET_NAME="east"
BASE_DIR=$(pwd) # bin dir
FSC_PATH="/shared/home/ynie/work2/software/fsc28_linux64"
EASY_SFS_SCRIPT_MODIFIED="./scripts/easySFS_withoutMultiSFS.py" # modified easySFS script to work without MultiSFS
PARALLEL_JOBS=50

# Set folders
DEMO_FILES="${BASE_DIR}/../output/demo_files"
EASY_SFS_DIR="${BASE_DIR}/../output/easySFS"
SUMMARY_DIR="${DEMO_FILES}/${SET_NAME}_summary"
BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_boot"
BLOCK_BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_block_boot"


myfunc_get_boot_blocks() {

    # Get the best run scenario name and run number
    SUMMARY_FILE="${SUMMARY_DIR}/Best_run_for_all.summary.txt"
    if [ ! -f "$SUMMARY_FILE" ]; then echo "Error: Summary file not found: $SUMMARY_FILE"; return 1; fi

    BEST_SCENARIO=$(head -n 2 "$SUMMARY_FILE" | tail -n 1 | cut -f1)
    BEST_RUN_NUM=$(head -n 2 "$SUMMARY_FILE" | tail -n 1 | cut -f2)
    echo "  -> Best Scenario: $BEST_SCENARIO | Run: $BEST_RUN_NUM" # e.g., eastG6D1 | run48

    # Create folders
    mkdir -p "$BOOT_DIR"
    mkdir -p "$BLOCK_BOOT_DIR/blocks"

    # Copy the best run files to the bootstrapping folder
    SRC_DIR="${DEMO_FILES}/${SET_NAME}"
    RUN_PATH=$(find "${SRC_DIR}" -type d -name "${BEST_SCENARIO}.${BEST_RUN_NUM}" | head -n 1)
    if [ -z "$RUN_PATH" ]; then echo "Error: Could not find run path for $BEST_SCENARIO.$BEST_RUN_NUM"; return 1; fi

    cp "${RUN_PATH}/${BEST_SCENARIO}/${BEST_SCENARIO}.pv" "$BOOT_DIR/"
    cp "${SRC_DIR}/Group_"*"/$BEST_SCENARIO.est" "$BOOT_DIR/"
    cp "${SRC_DIR}/Group_"*"/$BEST_SCENARIO.tpl" "$BOOT_DIR/"
    cp "${RUN_PATH}/"*_joint*.obs "$BOOT_DIR/"

    # Save bootstrap summary
    head -n1 "${RUN_PATH}/${BEST_SCENARIO}/${BEST_SCENARIO}.pv" > "$BOOT_DIR/$BEST_SCENARIO.bootstrap_summary.txt"

    # Seperate VCF header and body for block bootstrapping
    INPUT_VCF="${EASY_SFS_DIR}/${SET_NAME}.M5F1.syn_inter_LDprune.vcf"
 
    grep "^#" "$INPUT_VCF" > "$BLOCK_BOOT_DIR/$SET_NAME.vcf.header"
    grep -v "^#" "$INPUT_VCF" > "$BLOCK_BOOT_DIR/$SET_NAME.vcf.body"


    # Calculate number of sites per block (100 blocks total)
    siteCount=$(wc -l < "${BLOCK_BOOT_DIR}/${SET_NAME}.vcf.body")
    sitesPerBlock=$(( siteCount / 100 ))

    echo "  -> Total sites: $siteCount | Sites per block: $sitesPerBlock"

    split -l "$sitesPerBlock" -d "$BLOCK_BOOT_DIR/$SET_NAME.vcf.body" \
          "$BLOCK_BOOT_DIR/blocks/$SET_NAME.sites."
}


# Make bootstrapped SFS files
myfunc_bseasySFS() {
    BS_ID=$1 # $1 is the bootstrap number, e.g., 01, 02, ...

    # Make a new folder for each bootstrapping iteration:
    OUT_DIR="${BLOCK_BOOT_DIR}/${BS_ID}"
    mkdir -p "$OUT_DIR"

    TEMP_VCF="${OUT_DIR}/${SET_NAME}.bs.${BS_ID}.vcf"

    # Write header to new bootstrapped vcf file
    cat "${BLOCK_BOOT_DIR}/${SET_NAME}.vcf.header" > "$TEMP_VCF"

    # Randomly sampling and add 100 blocks
    echo "Generating bootstrapped VCF for BS_ID: $BS_ID"
    find "${BLOCK_BOOT_DIR}/blocks/" -name "${SET_NAME}.sites.*" | \
    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"${BS_ID}" -nosalt </dev/zero 2>/dev/null) -n 100 -r | \
    xargs cat >> "$TEMP_VCF"

    # Run easySFS
    echo "Running easySFS for bootstrapped VCF: $TEMP_VCF"
    
    # module load conda
    # . activate easySFS
    "$EASY_SFS_SCRIPT_MODIFIED" \
        -i "$TEMP_VCF" \
        -p "${EASY_SFS_DIR}/${SET_NAME}.group.set.txt" \
        --order betu,pash,SandCN,SandCNSW \
        -a -f \
        --proj 16,34,54,20 \
        -o "$OUT_DIR" \
        --prefix "$SET_NAME" > "${OUT_DIR}/easySFS.log" 2>&1
}
export -f myfunc_bseasySFS


## Run the functions
myfunc_get_boot_blocks
if [ $? -ne 0 ]; then echo "Step- Get boot blocks failed. Exiting."; exit 1; fi


bs=($(seq -w 1 50)) # Generate 50 shuffled vcf files
export SET_NAME
export BASE_DIR
export FSC_PATH
export EASY_SFS_SCRIPT_MODIFIED
export PARALLEL_JOBS
export DEMO_FILES
export EASY_SFS_DIR
export SUMMARY_DIR
export BOOT_DIR
export BLOCK_BOOT_DIR

bs=($(seq -w 1 50)) # Generate 50 shuffled vcf files
module load parallel
parallel -j 50 myfunc_bseasySFS ::: "${bs[@]}"

# export -f myfunc_bssbatach

# bs=($(seq -w 1 50)) # e.g., 01, 02, ..., 50

# for i in "${bs[@]}"; do
#     myfunc_bssbatch "$i"
# done


# mkdir -p "${BLOCK_BOOT_DIR}/summary" # Ensure summary directory exists
# perl ./scripts/fsc-selectbestrun.boot.pl \
#     "${BLOCK_BOOT_DIR}/*/*.run*/*/*.bestlhoods" \
#     "${BLOCK_BOOT_DIR}/summary"

# Rscript ./scripts/boot_intreval.R

# step_fsc_run