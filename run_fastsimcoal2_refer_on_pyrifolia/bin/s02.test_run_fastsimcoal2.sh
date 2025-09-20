#!/bin/bash

#SBATCH -J s02.test_run_fastsimcoal2.sh
#SBATCH -o s02.test_run_fastsimcoal2.sh.out
#SBATCH -e s02.test_run_fastsimcoal2.sh.err
#SBATCH -c 10
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=4
#SBATCH --mem=2G

# @File     :   s02.test_run_fastsimcoal2.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Time(Europe/Paris):   2025/09/19 16:06:36
#Description:
	# 

# step_fsc_test() { # fastsimcoal2 simulation, model test
#     mkdir -p ../output/demo_files/east
#     # a test run for east folder, the script will only run once with
#     # 100 simulation, and make the plot. It's better test them one by one,
#     # The running time is very short, so it won't take up time, but the
#     # output/error and the demo plots can help identify errors in tpl and est.
#     # I modified the ParFileViewer.r script to ParFileViewer.logname.r here
#     # to scale the y-axis in log10 and accept the pop names parameter.
#     # If your model fits better with the original version, please use the original script.
#     module load r/4.3.1
#     module load parallel
#     PATH=/shared/home/ynie/work2/software/fsc28_linux64/:$PATH
#     set="east"
#     workpath="/shared/home/ynie/work/pear/run_fastsimcoal2_refer_on_pyrifolia/bin"
#     Group="Group_G0"
#     prefix_list=("eastO1G0" "eastO2G0" "eastO3G0" "eastO4G0" \
#                  "eastT1G0" "eastT2G0" "eastT3G0" "eastT4G0" \
#                  "eastT5G0" "eastT6G0" "eastT7G0" "eastT8G0" )
	
#     for prefix in "${prefix_list[@]}"; do {
#         mkdir ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/easySFS/fastsimcoal2/*_joint*.obs \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.tpl \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.est \
#             ../output/demo_files/east/$Group/${prefix}.test_run

#         cd ../output/demo_files/east/$Group/${prefix}.test_run
#         for file in *.obs; do {
#             newname=$(echo "$file" | sed "s/$set/$prefix/g")
#             mv "$file" "$newname"
#         }; done
#         fsc28 -t ${prefix}.tpl -e ${prefix}.est \
#             -0 -C 10 -m \
#             -n 100 -L 40 -M -y 4 \
#             -q -c 50 -B 100
#         cd $prefix
#         Rscript $workpath/scripts/ParFileViewer.logname.r \
#             ${prefix}_maxL.par \
#             betu,pash,pyri_JP,Sand_CN >../../${prefix}_Viewer.log
#         cd $workpath

#     }; done

#     Group="Group_G1"
#         prefix_list=("eastO1G1" "eastO2G1" "eastO3G1" "eastO4G1" \
#                  "eastT1G1" "eastT2G1" "eastT3G1" "eastT4G1" \
#                  "eastT5G1" "eastT6G1" "eastT7G1" "eastT8G1" )
#     for prefix in "${prefix_list[@]}"; do {
#         mkdir ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/easySFS/fastsimcoal2/*_joint*.obs \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.tpl \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.est \
#             ../output/demo_files/east/$Group/${prefix}.test_run

#         cd ../output/demo_files/east/$Group/${prefix}.test_run
#         for file in *.obs; do {
#             newname=$(echo "$file" | sed "s/$set/$prefix/g")
#             mv "$file" "$newname"
#         }; done
#         fsc28 -t ${prefix}.tpl -e ${prefix}.est \
#             -0 -C 10 -m \
#             -n 100 -L 40 -M -y 4 \
#             -q -c 50 -B 100
#         cd $prefix
#         Rscript $workpath/scripts/ParFileViewer.logname.r \
#             ${prefix}_maxL.par \
#             betu,pash,pyri_JP,Sand_CN >../../${prefix}_Viewer.log
#         cd $workpath

#     }; done

#     Group="Group_G2"
# 	prefix_list=("eastT5G2" "eastT6G2" "eastT7G2" "eastT8G2" )
#     for prefix in "${prefix_list[@]}"; do {
#         mkdir ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/easySFS/fastsimcoal2/*_joint*.obs \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.tpl \
#             ../output/demo_files/east/$Group/${prefix}.test_run
#         cp ../output/demo_files/east/$Group/${prefix}.est \
#             ../output/demo_files/east/$Group/${prefix}.test_run

#         cd ../output/demo_files/east/$Group/${prefix}.test_run
#         for file in *.obs; do {
#             newname=$(echo "$file" | sed "s/$set/$prefix/g")
#             mv "$file" "$newname"
#         }; done
#         fsc28 -t ${prefix}.tpl -e ${prefix}.est \
#             -0 -C 10 -m \
#             -n 100 -L 40 -M -y 4 \
#             -q -c 50 -B 100
#         cd $prefix
#         Rscript $workpath/scripts/ParFileViewer.logname.r \
#             ${prefix}_maxL.par \
#             betu,pash,pyri_JP,Sand_CN >../../${prefix}_Viewer.log
#         cd $workpath
#     }; done

# }

# step_fsc_test


run_fsc_test() {
    Group="$1"
    prefix="$2"
    set="east"
    workpath="/shared/home/ynie/work/pear/run_fastsimcoal2_refer_on_pyrifolia/bin"
    PATH=/shared/home/ynie/work2/software/fsc28_linux64/:$PATH
    module load r/4.3.1
    module load parallel
    cd $workpath
    mkdir -p ../output/demo_files/east/$Group/${prefix}.test_run
    cp ../output/easySFS/fastsimcoal2/*_joint*.obs \
        ../output/demo_files/east/$Group/${prefix}.test_run
    cp ../output/demo_files/east/$Group/${prefix}.tpl \
        ../output/demo_files/east/$Group/${prefix}.test_run
    cp ../output/demo_files/east/$Group/${prefix}.est \
        ../output/demo_files/east/$Group/${prefix}.test_run

    cd ../output/demo_files/east/$Group/${prefix}.test_run
    for file in *.obs; do
        newname=$(echo "$file" | sed "s/${set}_/${prefix}_/g")
        mv "$file" "$newname"
    done
    fsc28 -t ${prefix}.tpl -e ${prefix}.est \
        -0 -C 10 -m \
        -n 100 -L 40 -M -y 4 \
        -q -c 50 -B 100
    cd $prefix
    Rscript $workpath/scripts/ParFileViewer.logname.r \
        ${prefix}_maxL.par \
        betu,pash,pyri_JP,Sand_CN >../../${prefix}_Viewer.log
    cd $workpath
}
export -f run_fsc_test

step_fsc_test() {
    workpath="/shared/home/ynie/work/pear/run_fastsimcoal2_refer_on_pyrifolia/bin"
    set="east"
    module load parallel
    parallel --jobs 4 run_fsc_test Group_G0 ::: \
        eastO1G0 eastO2G0 eastO3G0 eastO4G0 \
        eastT1G0 eastT2G0 eastT3G0 eastT4G0 \
        eastT5G0 eastT6G0 eastT7G0 eastT8G0 

    parallel --jobs 4 run_fsc_test Group_G1 ::: \
        eastO1G1 eastO2G1 eastO3G1 eastO4G1 \
        eastT1G1 eastT2G1 eastT3G1 eastT4G1 \
        eastT5G1 eastT6G1 eastT7G1 eastT8G1 

    parallel --jobs 4 run_fsc_test Group_G2 ::: \
        eastT5G2 eastT6G2 eastT7G2 eastT8G2

    wait
}

mkdir -p ../output/demo_files/east
cp -r ../input/demography_models_config/* ../output/demo_files/east

step_fsc_test

cd ../output/demo_files/east
mkdir -p PDFs.test_run
mv Group_*/*.test_run/*/*.pdf  PDFs.test_run
