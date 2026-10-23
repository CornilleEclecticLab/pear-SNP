#!/bin/bash
# By Xilong CHEN
# Create date: 2025-08-17
# Contact: chen_xilong@outlook.com

# 2025-08-21 updated for comm ref
# 2025-09-20 Yuqi updated for ppy ref
# 2026-09-01 Yuqi updated for IFB

mkdir -p ../output/

# four pop name
# betu
# pash
# pyri_JP
# Sand_CN

step_fsc_run() { # fastsimcoal2 simulation, run 50 times for each model 100000 simulations
    set="east"
    # workpath="/home/xchen/work/PearPop_Project_Yuqi"
    # workpath="/scratch/yn2515/pear_snp_tip/run_fastsimcoal2_refer_on_pyrifolia_4pop_b/bin"
    workpath=$(pwd)
    # Group_list=("Group_G0" "Group_G1" "Group_G2" "Group_G3")

    G_list=("G0" "G1")
    total_runs=50
    runs_per_job=5
    for g in "${G_list[@]}"; do
        prefix_list=()
        Group="Group_${g}"
        for d in {1..8}; do
            for i in "S"; do
                prefix="${set}${i}${d}${g}"
                if [ -f ../output/demo_files/$set/$Group/${prefix}.tpl ]; 
                then
                    prefix_list+=($prefix)
                fi
            done
        done
        echo $Group
        echo "${prefix_list[@]}"

        for prefix in "${prefix_list[@]}"; do
            for i in $(seq -f "%02g" 1 $total_runs); do
                echo "Prepareing $prefix run$i files"
                pwd
                rundir="../output/demo_files/${set}/$Group/${prefix}.run$i"
                if [ -d $rundir ]; then
                    echo "Error: Directory $rundir already exists."
                    exit 1
                    # rm -rf ../output/demo_files/${set}/$Group/${prefix}.run$i
                fi
                mkdir "$rundir"
                cp ../output/easySFS/fastsimcoal2/*_joint*.obs "$rundir"
                cp "../output/demo_files/$set/$Group/${prefix}.tpl" "$rundir"
                cp "../output/demo_files/$set/$Group/${prefix}.est" "$rundir"

                (
                cd "$rundir" || exit 1
                for file in *.obs; do
                    newname=$(echo "$file" | sed "s/${set}_/${prefix}_/g")
                    mv "$file" "$newname"
                done
                )
            done
        done
        # submit jobs
        for prefix in "${prefix_list[@]}"; do
        batch_id=0
        for batch_start in $(seq 1 $runs_per_job $total_runs); do
            batch_id=$((batch_id + 1))
            batch_end=$((batch_start + runs_per_job - 1))
                sbatch_dir="../output/demo_files/${set}/$Group"
                sbatch_file="${sbatch_dir}/${prefix}.batch$(printf '%02d' $batch_id).sbatch"



            # for i in $(seq -w 1 50); do
            #     cd ../output/demo_files/$set/$Group/${prefix}.run$i
                cat <<END >"$sbatch_file"
#!/bin/bash
#SBATCH -J ${prefix}.batch$(printf '%02d' $batch_id)
#SBATCH -o ${prefix}.batch$(printf '%02d' $batch_id).%J.out
#SBATCH -e ${prefix}.batch$(printf '%02d' $batch_id).%J.err

PATH=/shared/home/ynie/work2/software/fsc28_linux64:\$PATH

for i in \$(seq -f "%02g" $batch_start $batch_end); do
    (
        rundir="${workpath}/../output/demo_files/${set}/$Group/${prefix}.run\$i"
    cd "\$rundir" || exit 1
    fsc28 -t ${prefix}.tpl -e ${prefix}.est \\
        -0 -C 10 -m \\
        -n 100000 -L 40 -M -y 4 \\
        -q -c 10 -B 20
    )
done
END
                sbatch -c 10 "$sbatch_file"
            done
        done
    done

}

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

step_fsc_run