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
FSC_PATH="/shared/projects/pear_snp2/software/fsc28_linux64"
EASY_SFS_SCRIPT_MODIFIED="easySFS_withoutMultiSFS.py" # modified easySFS script to work without MultiSFS
PARALLEL_JOBS=50

# Set folders
DEMO_FILES="${BASE_DIR}/../output/demo_files"
EASY_SFS_DIR="${BASE_DIR}/../output/easySFS"
SUMMARY_DIR="${DEMO_FILES}/${SET_NAME}_summary"
BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_boot"
BLOCK_BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_block_boot"      # old results
INIT_BOOT_DIR="${DEMO_FILES}/${SET_NAME}_Bestrun_block_boot_init"  # rerun: --initValues


# # Run fastsimcoal2 on bootstrapped SFS files
myfunc_bssbatch() {
    BS_ID=$1

    # Get the best run scenario name
    SUMMARY_FILE="${SUMMARY_DIR}/Best_run_for_all.summary.txt"
    BEST_SCENARIO=$(head -n 2 "$SUMMARY_FILE" | tail -n 1 | cut -f1)

    # Parameters
    RE_COUNT=5
    RUNS_PER_JOB=1
    SIM_COUNT=100000
    EMC_CYCLE=40
    CPU_PER_JOB=10
    HOURS_PER_RUN=4

    SRC_SUB_DIR="${BLOCK_BOOT_DIR}/${BS_ID}"
    BS_SUB_DIR="${INIT_BOOT_DIR}/${BS_ID}"
    OBS_SOURCE_DIR="${BS_SUB_DIR}/fastsimcoal2"

    if [ ! -d "${SRC_SUB_DIR}/fastsimcoal2" ]; then
        echo "Error: expected easySFS output not found: ${SRC_SUB_DIR}/fastsimcoal2 (bootstrap ${BS_ID})" >&2
    exit 1
    fi
    if [ ! -f "${BOOT_DIR}/${BEST_SCENARIO}.pv" ]; then
        echo "Error: best-run .pv not found: ${BOOT_DIR}/${BEST_SCENARIO}.pv" >&2
    exit 1
    fi
    if compgen -G "${BS_SUB_DIR}/${SET_NAME}.run*" > /dev/null; then
        echo "Error: ${BS_SUB_DIR} already has run dirs; refusing to overwrite (bootstrap ${BS_ID})" >&2
    exit 1
    fi

    # Copy the bootstrap SFS (easySFS output) from the old block_boot dir
    mkdir -p "$BS_SUB_DIR"
    cp -rp "${SRC_SUB_DIR}/fastsimcoal2" "$BS_SUB_DIR/"
    cp -p "${SRC_SUB_DIR}/datadict.txt" "${SRC_SUB_DIR}/easySFS.log" "$BS_SUB_DIR/" 2>/dev/null || true

    for i in $(seq -f "%02g" 1 ${RE_COUNT}); do
        RUN_DIR="${BS_SUB_DIR}/${SET_NAME}.run${i}"
        mkdir -p "$RUN_DIR"

        # Copy the bootstrapped .obs files
        cp "${OBS_SOURCE_DIR}/"*_joint*.obs "$RUN_DIR/"
        # Copy the best scenario tpl and est files
        cp "${BOOT_DIR}/${BEST_SCENARIO}.tpl" "$RUN_DIR/${SET_NAME}.tpl"
        cp "${BOOT_DIR}/${BEST_SCENARIO}.est" "$RUN_DIR/${SET_NAME}.est"
        # Start from the ML estimates of the original data
        cp "${BOOT_DIR}/${BEST_SCENARIO}.pv" "$RUN_DIR/${SET_NAME}.pv"
    done

    batch_id=0
    for batch_start in $(seq 1 ${RUNS_PER_JOB} ${RE_COUNT}); do
        batch_id=$((batch_id + 1))
        batch_end=$((batch_start + RUNS_PER_JOB - 1))
        if [ "$batch_end" -gt "$RE_COUNT" ]; then batch_end=$RE_COUNT; fi
        n_runs_this_batch=$((batch_end - batch_start + 1))

        SCRIPT_NAME="${BS_SUB_DIR}/bs${BS_ID}.batch$(printf '%02d' $batch_id).sbatch"
        cat <<END > "$SCRIPT_NAME"
#!/bin/bash
#SBATCH -J bs${BS_ID}.b$(printf '%02d' $batch_id)
#SBATCH -o bs${BS_ID}.b$(printf '%02d' $batch_id).%J.out
#SBATCH -e bs${BS_ID}.b$(printf '%02d' $batch_id).%J.err
#SBATCH -D ${BS_SUB_DIR}
#SBATCH -n 1
#SBATCH -c ${CPU_PER_JOB}
#SBATCH --mem=4G
#SBATCH --time=${HOURS_PER_RUN}:00:00

PATH=${FSC_PATH}:\$PATH

for i in \$(seq -f "%02g" ${batch_start} ${batch_end}); do
    (
        rundir="${BS_SUB_DIR}/${SET_NAME}.run\$i"
        cd "\$rundir" || exit 1
        fsc28 -t ${SET_NAME}.tpl -e ${SET_NAME}.est \\
        -0 -C 10 -m \\
        -n ${SIM_COUNT} -L ${EMC_CYCLE} -M -y 4 \\
        --initValues ${SET_NAME}.pv \\
        -q -c ${CPU_PER_JOB} -B 20
    )
done
END
        # Submit the job
        sbatch -c ${CPU_PER_JOB} "$SCRIPT_NAME"
    done
}

export -f myfunc_bssbatch


### Run the functions
bs=(${BS_LIST:-$(seq -w 1 50)}) # e.g., 01, 02, ..., 50; test one: BS_LIST=01 bash s06.IFB.bootstrap.sh

for i in "${bs[@]}"; do
    myfunc_bssbatch "$i"
done


# mkdir -p "${BLOCK_BOOT_DIR}/summary" # Ensure summary directory exists
# perl ./scripts/fsc-selectbestrun.boot.pl \
#     "${BLOCK_BOOT_DIR}/*/*.run*/*/*.bestlhoods" \
#     "${BLOCK_BOOT_DIR}/summary"

# Rscript ./scripts/boot_intreval.R

# step_fsc_run