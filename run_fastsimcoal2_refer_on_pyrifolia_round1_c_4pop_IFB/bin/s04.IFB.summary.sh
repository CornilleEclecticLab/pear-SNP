#!/bin/bash
# By Xilong CHEN
# Create date: 2025-08-17
# Contact: chen_xilong@outlook.com

# 2025-08-21 updated for comm ref
# 2025-09-20 Yuqi updated for ppy ref
# 2026-09-01 Yuqi updated for IFB

mkdir -p ../output/


step_summary() {
    set="east"
    workpath=$(pwd)
    mkdir -p ../output/demo_files/${set}_summary
    perl ./scripts/fsc-selectbestrun.pl \
        "../output/demo_files/${set}/*/*.run*/*/*.bestlhoods" \
        ../output/demo_files/${set}_summary
    cat ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt

    # Scenario        run_nu  MaxEstLhood     AIC     deltaL  best_of_what
    # eastT6G2        run23   -1716002.013    7902507.30936315        101617.161      best_of_MaxEstLhood
    # eastT6G2        run23   -1716002.013    7902507.30936315        101617.161      best_of_AIC

    BEST_SCENARIO=$(head -n 2 ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f1)
    echo "Best scenario is: $BEST_SCENARIO"  # e.g., eastT6G2
    BEST_RUN_NUM=$(head -n 2 ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f2)
    echo "Best run number is: $BEST_RUN_NUM"  # e.g., run23

    mkdir -p ../output/demo_files/${set}_Bestrun_bestlhoods_file
    cp ../output/demo_files/${set}/Group_G*/${BEST_SCENARIO}.${BEST_RUN_NUM}/${BEST_SCENARIO}/*.bestlhoods \
        ../output/demo_files/${set}_Bestrun_bestlhoods_file/
    cp ../output/demo_files/${set}/Group_G*/${BEST_SCENARIO}.${BEST_RUN_NUM}/${BEST_SCENARIO}/*.par \
        ../output/demo_files/${set}_Bestrun_bestlhoods_file/

    cd ../output/demo_files/${set}_Bestrun_bestlhoods_file/

    module load r/4.3.1
    Rscript $workpath/scripts/ParFileViewer.logname.r \
        ${BEST_SCENARIO}_maxL.par \
        betu,pash,Sand_CN,Sand_CN-SW # cauc,comD,comP,pyra
    cd $workpath

# compare the scenarios
    Rscript ./scripts/fsc_MaxEsthood_boxplot.R
    # Rscript ./scripts/fsc_best_confidence_interval.R

}

# step_fsc_run
step_summary
