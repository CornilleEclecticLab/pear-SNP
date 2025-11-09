#!/bin/bash
# By Xilong CHEN
# Create date: 2025-08-17
# Contact: chen_xilong@outlook.com

# 2025-09-20 Yuqi updated for ppy ref

mkdir -p ../output/


step_fsc_run() { # fastsimcoal2 simulation, run 50 times for each model 100000 simulations
    set="east"
    # use caller's workpath if exported, otherwise fall back to current directory
    : "${workpath:=$(pwd)}"

    # Group_list=("Group_G0" "Group_G1" "Group_G2" "Group_G3")
    G_list=("G0" "G1" "G2")
    for g in "${G_list[@]}"; do
        prefix_list=()
        Group="Group_${g}"
        for d in {1..4}; do
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
            for i in $(seq -w 1 50); do
                echo "Prepareing $prefix run$i files"
                pwd
                if [ -d ../output/demo_files/${set}/$Group/${prefix}.run$i ]; then
                    echo "Error: Directory ../output/demo_files/${set}/$Group/${prefix}.run$i already exists."
                    exit 1
                    # rm -rf ../output/demo_files/${set}/$Group/${prefix}.run$i
                fi
                mkdir ../output/demo_files/${set}/$Group/${prefix}.run$i
                cp ../output/easySFS/fastsimcoal2/*_joint*.obs \
                    ../output/demo_files/$set/$Group/${prefix}.run$i
                cp ../output/demo_files/$set/$Group/${prefix}.tpl \
                    ../output/demo_files/$set/$Group/${prefix}.run$i
                cp ../output/demo_files/$set/$Group/${prefix}.est \
                    ../output/demo_files/$set/$Group/${prefix}.run$i

                cd ../output/demo_files/$set/$Group/${prefix}.run$i
                for file in *.obs; do
                    newname=$(echo "$file" | sed "s/${set}_/${prefix}_/g")
                    mv "$file" "$newname"
                done
                cd $workpath || exit
            done
        done
        # submit jobs
        for prefix in "${prefix_list[@]}"; do
            for i in $(seq -w 1 50); do
                cd ../output/demo_files/$set/$Group/${prefix}.run$i
                cat <<END >$prefix.run$i.sbatch
#!/bin/bash
#SBATCH -J $prefix.run$i
#SBATCH -o $prefix.run$i.%J.out
#SBATCH -e $prefix.run$i.%J.err

. $workpath/load_fastsimcoal.sh

fsc28 -t ${prefix}.tpl -e ${prefix}.est \
    -0 -C 10 -m \
    -n 100000 -L 40 -M -y 4 \
    -q -c 10 -B 20
END
                sbatch -c 10 ${prefix}.run$i.sbatch
                cd "$workpath"
            done
        done
    done

}

step_summary() {
    set="east"
    : "${workpath:=$(pwd)}"

    mkdir ../output/demo_files/${set}_summary

    # This perl script will extract the best likelihoods from all runs and summarize them
    perl ./scripts/fsc-selectbestrun.pl \
        "../output/demo_files/${set}/*/*.run*/*/*.bestlhoods" \
        ../output/demo_files/${set}_summary
    cat ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt

# It will print something like this:
# Scenario        run_nu  MaxEstLhood     AIC     deltaL  best_of_what
# eastS2G1        run23   -2132725.395    9821585.40395368        125007.2        best_of_MaxEstLhood
# eastS2G1        run23   -2132725.395    9821585.40395368        125007.2        best_of_AIC
    
    # Capture the best scenario and best run number
    BST_SCENARIO=$(head -n 2 ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f1)
    echo "Best scenario is: $BEST_SCENARIO"
    BEST_RUN_NUM=$(head -n 2 ../output/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f2)

    mkdir ../output/demo_files/${set}_Bestrun_bestlhoods_file
    cp ../output/demo_files/${set}/Group_G*/${BEST_SCENARIO}.${BEST_RUN_NUM}/${BEST_SCENARIO}/*.bestlhoods \
        ../output/demo_files/${set}_Bestrun_bestlhoods_file/
    cp ../output/demo_files/${set}/Group_G*/${BEST_SCENARIO}.${BEST_RUN_NUM}/${BEST_SCENARIO}/*.par \
        ../output/demo_files/${set}_Bestrun_bestlhoods_file/

    cd ../output/demo_files/${set}_Bestrun_bestlhoods_file/

    . s00.load_r.sh

    # Visualize the parameter estimates
    Rscript $workpath/scripts/ParFileViewer.logname.r \
        ${set}S2G1_maxL.par \
        betu,pash,Sand_CN,Sand_CN-SW
    cd $workpath

    # Visualize the AIC values, the plot will be saved in the folder ../output/demo_files/${set}_summary
    Rscript ./scripts/fsc_MaxEsthood_boxplot.commbined.R


}

step_fsc_run
