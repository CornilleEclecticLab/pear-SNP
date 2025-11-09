#!/bin/bash

#SBATCH -J s02.test_run_fastsimcoal2.sh
#SBATCH -o s02.test_run_fastsimcoal2.sh.%J.out
#SBATCH -e s02.test_run_fastsimcoal2.sh.%J.err
#SBATCH -c 8
#SBATCH --mem=8G

# @File     :   s02.test_run_fastsimcoal2.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Time(Europe/Paris):   2025/09/19 16:06:36
#Description:
	# This script is used to run fastsimcoal2 test runs for demo files, 
    # generate PDF files for each test run,
    # and copy .obs, .tpl and .est files to proper position for simulation runs, i.e., s03.Jubail.run_fastsimcoal2.sh



run_fsc_test() {
    Group="$1"
    prefix="$2"
    set="east"
    # use caller's workpath if exported, otherwise fall back to current directory
    : "${workpath:=$(pwd)}"  # bin directory
    . load_fastsimcoal2.sh
    . load_r.sh
    cd "$workpath"
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
        betu,pash,Sand_CN,Sand_CN-SW >../../${prefix}_Viewer.log
    cd $workpath
}
export -f run_fsc_test

step_fsc_test() {
    : "${workpath:=$(pwd)}"
    set="east"
    . s00.load_parallel.sh

    parallel --jobs 4 run_fsc_test Group_G0 ::: \
        eastS1G0 eastS2G0 eastS3G0 eastS4G0

    parallel --jobs 4 run_fsc_test Group_G1 ::: \
        eastS1G1 eastS2G1 eastS3G1 eastS4G1

    wait
}


mkdir -p ../output/demo_files/east
cp -r ../input/demography_models_config/* ../output/demo_files/east

# Launch fastsimcoal2 test runs
step_fsc_test

# Collect all PDF files to one folder "../output/demo_files/east/PDFs.test_run"
cd ../output/demo_files/east
mkdir -p PDFs.test_run
mv Group_*/*.test_run/*/*.pdf  PDFs.test_run
