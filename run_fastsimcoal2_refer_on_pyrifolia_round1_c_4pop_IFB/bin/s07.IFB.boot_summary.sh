#!/bin/bash

# @File     :   s04.Jubail.bootstrap.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi modified from Xilong's scripts
# @Time(CET):   2025/11/20 09:36:34
#Description:
    # 


set -euo pipefail

### Configure
SET_NAME="east"
BASE_DIR=$(pwd) # bin dir
EASY_SFS_SCRIPT_MODIFIED="easySFS_withoutMultiSFS.py" # modified easySFS script to work without MultiSFS
PARALLEL_JOBS=50

# Set folders
DEMO_FILES="${BASE_DIR}/../output/demo_files"
EASY_SFS_DIR="${BASE_DIR}/../output/easySFS"
SUMMARY_DIR="${DEMO_FILES}/${SET_NAME}_summary"
BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_boot"
BLOCK_BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_block_boot_init"


mkdir -p "${BLOCK_BOOT_DIR}/summary" # Ensure summary directory exists
perl ./scripts/fsc-selectbestrun.boot.pl \
    "${BLOCK_BOOT_DIR}/*/*.run*/*/*.bestlhoods" \
    "${BLOCK_BOOT_DIR}/summary"
