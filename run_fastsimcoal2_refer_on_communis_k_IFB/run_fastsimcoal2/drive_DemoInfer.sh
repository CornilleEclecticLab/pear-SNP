#!/bin/bash
# By Xilong CHEN
# Create date: 2025-08-17
# Contact: chen_xilong@outlook.com

# 2025-08-21 updated for comm ref

# 2026-08-10 update for running on IFB

mkdir -p ./analysis/DemoInfer/

# four pop name
# cauc
# comD
# comP
# pyra

# sort -k1,1 -k4,4n ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.gff | grep -v "^$" |
#    awk '{if($5 < $4) $5=$4; print}' OFS="\t" \
#        >./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff

# bgzip ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff
# tabix -p gff ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff.gz

# zcat ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.sorted.gff.gz | sed 's/;;/;/g' | bgzip \
#     >./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff.gz
# tabix -p gff ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff.gz

step_PrepareVcf() {
    mkdir -p ./analysis/DemoInfer/PrepareVcf

    bcftools filter -e 'F_MISSING > 0.2' \
        ./data/fsc_input_250816/pear_Jul2024_ref_comm.Europe_noAdmix.Combine_chr.maf005.vcf.gz -Oz \
        -o ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz
    tabix -p vcf ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz

    # VEP
    mkdir -p ./analysis/DemoInfer/PrepareVcf/syn_inter
    vep \
        --input_file ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz \
        --fasta ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fasta \
        --gff ./data/fsc_input_250816/PyrusCommunis_BartlettDHv2.0.fixed.gff.gz \
        --species WestPear \
        --distance 2000 \
        --force_overwrite \
        --fork 8 \
        --output_file ./analysis/DemoInfer/PrepareVcf/syn_inter/west.M5F1.gz.f2k.annotated.vep.txt

    cat ./analysis/DemoInfer/PrepareVcf/syn_inter/west.M5F1.gz.f2k.annotated.vep.txt |
        grep -E 'synonymous_variant|intergenic_variant' |
        awk '{print $2}' |
        sed 's/:/\t/g' \
            >./analysis/DemoInfer/PrepareVcf/syn_inter/west.M5F1.gz.f2k.syn_inter.vep.txt

    bcftools view -R ./analysis/DemoInfer/PrepareVcf/syn_inter/west.M5F1.gz.f2k.syn_inter.vep.txt \
        ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz -Oz \
        -o ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter.vcf.gz
    tabix -p vcf ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter.vcf.gz

    # ld

    mkdir -p ./analysis/DemoInfer/PrepareVcf/plink_LDprune
    plink --vcf ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter.vcf.gz \
        --make-bed \
        --indep-pairwise 50 5 0.2 \
        --const-fid \
        --out ./analysis/DemoInfer/PrepareVcf/plink_LDprune/west.M5F1.syn_inter \
        --set-missing-var-ids \@\:# \
        --keep-allele-order \
        --allow-extra-chr

    # prune.in to site
    cat ./analysis/DemoInfer/PrepareVcf/plink_LDprune/west.M5F1.syn_inter.prune.in |
        sed 's/:/\t/g' \
            >./analysis/DemoInfer/PrepareVcf/plink_LDprune/west.M5F1.syn_inter.prune.in.site

    # keep prune.in site
    bcftools view \
        -R ./analysis/DemoInfer/PrepareVcf/plink_LDprune/west.M5F1.syn_inter.prune.in.site \
        ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter.vcf.gz \
        -Oz -o ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter_LDprune.vcf.gz
    tabix -p vcf ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter_LDprune.vcf.gz
    # 135677
}

step_easySFS() {
    module load conda
    . activate easySFS
    module load r/4.3.1
    module load python

    mkdir -p ./analysis/DemoInfer/easySFS
    # prepare a pure vcf file

    gunzip -c ./analysis/DemoInfer/PrepareVcf/west.M5F1.syn_inter_LDprune.vcf.gz \
        >./analysis/DemoInfer/easySFS/west.M5F1.syn_inter_LDprune.vcf
    # prepare a pop/group file
    sed -e 's/comm_Dessert/comD/g' -e 's/comm_Perry/comP/g' \
        ./data/fsc_input_250816/s01.sample_tab_population.western.txt \
        >./analysis/DemoInfer/easySFS/west.group.set.txt

    # easySFS part, here the vcf is not the RADvcf, so use -a to consider the full sites
    # Preview
    # jubail py3 with pandas, scipy, Aphid env

    python3 ./bin/easySFS/easySFS.py \
        -i ./analysis/DemoInfer/easySFS/west.M5F1.syn_inter_LDprune.vcf \
        -p ./analysis/DemoInfer/easySFS/west.group.set.txt \
        --order cauc,comD,comP,pyra \
        -a \
        --preview \
        >./analysis/DemoInfer/easySFS/preview.output.txt
    # find the candidate project from output of preview
    ## convert output for R plot
    perl ./scripts/easySFSpreview_output_format.pl \
        ./analysis/DemoInfer/easySFS/preview.output.txt \
        >./analysis/DemoInfer/easySFS/preview.output.Rdata.txt
    ## plot them and output the max value
    Rscript ./scripts/easySFSpreview_plot.R \
        >./analysis/DemoInfer/easySFS/preview.output.proj.result.txt
    ## the result is
    # cauc,comD,comP,pyra
    # 42,28,22,18
    # Here the projection should also be used for tpl and est file
    # easySFS running

    ## the multiSFS take too much time and memory for big pop size, here I
    ## remove this function
    ## order cauc,comD,comP,pyra
    ## 42,28,22,18
    python3 ./bin/easySFS/easySFS_withoutM.py \
        -i ./analysis/DemoInfer/easySFS/west.M5F1.syn_inter_LDprune.vcf \
        -p ./analysis/DemoInfer/easySFS/west.group.set.txt \
        --order cauc,comD,comP,pyra \
        -a -f \
        --proj 42,28,22,18 \
        -o ./analysis/DemoInfer/easySFS \
        --prefix west

}

step_fsc_test() { # fastsimcoal2 simulation, model test
    mkdir -p ./analysis/DemoInfer/demo_files
    mkdir -p ./analysis/DemoInfer/demo_files/west
    # a test run for west folder, the script will only run once with
    # 100 simulation, and make the plot. It's better test them one by one,
    # The running time is very short, so it won't take up time, but the
    # output/error and the demo plots can help identify errors in tpl and est.
    # I modified the ParFileViewer.r script to ParFileViewer.logname.r here
    # to scale the y-axis in log10 and accept the pop names parameter.
    # If your model fits better with the original version, please use the original script.

    module load r/4.3.1
    PATH=/shared/home/ynie/work2/software/fsc28_linux64/:$PATH

    set="west"
    workpath=$(pwd)
    Group="Group_G0"
    prefix_list=("westG0D1" "westG0D2" "westG0D3" "westG0D4" "westG0D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G1"
    prefix_list=("westG1D1" "westG1D2" "westG1D3" "westG1D4" "westG1D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G2"
    prefix_list=("westG2D1" "westG2D2" "westG2D3" "westG2D4" "westG2D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G3"
    prefix_list=("westG3D1" "westG3D2" "westG3D3" "westG3D4" "westG3D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G4"
    prefix_list=("westG4D1" "westG4D2" "westG4D3" "westG4D4" "westG4D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G5"
    prefix_list=("westG5D1" "westG5D2" "westG5D3" "westG5D4" "westG5D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G6"
    prefix_list=("westG6D1" "westG6D2" "westG6D3" "westG6D4" "westG6D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G7"
    prefix_list=("westG7D1" "westG7D2" "westG7D3" "westG7D4" "westG7D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G8"
    prefix_list=("westG8D1" "westG8D2" "westG8D3" "westG8D4" "westG8D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G9"
    prefix_list=("westG9D1" "westG9D2" "westG9D3" "westG9D4" "westG9D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G10"
    prefix_list=("westG10D1" "westG10D2" "westG10D3" "westG10D4" "westG10D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done

    Group="Group_G11"
    prefix_list=("westG11D1" "westG11D2" "westG11D3" "westG11D4" "westG11D5")
    for prefix in "${prefix_list[@]}"; do {
        mkdir -p ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
            ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run

        cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
        for file in *.obs; do {
            newname=$(echo "$file" | sed "s/$set/$prefix/g")
            mv "$file" "$newname"
        }; done
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \
            -0 -C 10 -m \
            -n 100 -L 40 -M -y 4 \
            -q -c 50 -B 100
        cd $prefix
        Rscript $workpath/scripts/ParFileViewer.logname.r \
            ${prefix}_maxL.par \
            cauc,comD,comP,pyra >../../${prefix}_Viewer.log
        cd $workpath

    }; done


}

step_fsc_run() { # fastsimcoal2 simulation, run 50 times for each model 100000 simulations
    set="west"
    workpath=$(pwd)
    # Group_list=("Group_G0" "Group_G1" "Group_G2" "Group_G3")

    G_list=("G0" "G1" "G2" "G3" "G4" "G5" "G6" "G7" "G8" "G9" "G10" "G11")
    # G_list=("G10" "G11")
    total_runs=50
    runs_per_job=5

    for g in "${G_list[@]}"; do
        prefix_list=()
        Group="Group_${g}"
        for d in 1 2 3 4; do
        #for d in 5; do
            prefix_list+=("west${g}D${d}")
        done
        echo $Group
        echo "${prefix_list[@]}"

        for prefix in "${prefix_list[@]}"; do
            for i in $(seq -f "%02g" 1 $total_runs); do
                rundir="./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i"
                mkdir -p "$rundir"
                cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs "$rundir"
                cp "./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl" "$rundir"
                cp "./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est" "$rundir"
                (
                    cd "$rundir" || exit 1
                    for file in *.obs; do
                        newname=$(echo "$file" | sed "s/$set/$prefix/g")
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
                sbatch_dir="./analysis/DemoInfer/demo_files/west/$Group"
                sbatch_file="${sbatch_dir}/${prefix}.batch$(printf '%02d' $batch_id).sbatch"

                cat <<END >"$sbatch_file"
#!/bin/bash
#SBATCH -J ${prefix}.batch$(printf '%02d' $batch_id)
#SBATCH -o ${prefix}.batch$(printf '%02d' $batch_id).%J.out
#SBATCH -e ${prefix}.batch$(printf '%02d' $batch_id).%J.err


PATH=/shared/home/ynie/work2/software/fsc28_linux64:\$PATH

for i in \$(seq -f "%02g" $batch_start $batch_end); do
        rundir="${workpath}/analysis/DemoInfer/demo_files/west/$Group/${prefix}.run\$i"
        cd "\$rundir" || exit 1
        fsc28 -t ${prefix}.tpl -e ${prefix}.est \\
            -0 -C 10 -m \\
            -n 100000 -L 40 -M -y 4 \\
            -q -c 10 -B 20
done
END
                sbatch -c 10 "$sbatch_file"
            done
        done
    done

}

step_summary() {
    set="west"
    module load r/4.3.1

    mkdir -p ./analysis/DemoInfer/demo_files/west_summary
    perl ./scripts/fsc-selectbestrun.pl \
        "./analysis/DemoInfer/demo_files/west/*/*.run*/*/*.bestlhoods" \
        ./analysis/DemoInfer/demo_files/west_summary
    cat ./analysis/DemoInfer/demo_files/west_summary/Best_run_for_all.summary.txt
    # Scenario        run_nu  MaxEstLhood     AIC     deltaL  best_of_what
    # westG6D1        run48   -1833886.714    8445380.41979247        96438.9639999999        best_of_MaxEstLhood
    # westG6D1        run48   -1833886.714    8445380.41979247        96438.9639999999        best_of_AIC

    # Get the best run scenario name and run number
    BEST_SCENARIO=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f1)
    echo "Best scenario is: $BEST_SCENARIO"  # e.g., westG6D1

    BEST_RUN_NUM=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f2)
    echo "Best run number is: $BEST_RUN_NUM"  # e.g., run48
    mkdir -p ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file

    cp ./analysis/DemoInfer/demo_files/west/Group_G*/$BEST_SCENARIO.$BEST_RUN_NUM/$BEST_SCENARIO/*.bestlhoods \
        ./analysis/DemoInfer/demo_files/${set}_Bestrun_bestlhoods_file/
    cp ./analysis/DemoInfer/demo_files/west/Group_G*/$BEST_SCENARIO.$BEST_RUN_NUM/$BEST_SCENARIO/*.par \
        ./analysis/DemoInfer/demo_files/${set}_Bestrun_bestlhoods_file/

    # workpath="/home/xchen/work/PearPop_Project_Yuqi"
    workpath=$(pwd)
    cd ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file/
    Rscript $workpath/scripts/ParFileViewer.logname.r \
        ${BEST_SCENARIO}_maxL.par \
        cauc,comD,comP,pyra
    cd $workpath

# compare the scenarios
    Rscript ./scripts/fsc_MaxEsthood_boxplot.R
    # Rscript ./scripts/fsc_best_confidence_interval.R

}

myfunc_get_boot_blocks() {
set="west"
# Get the best run scenario name and run number
BEST_SCENARIO=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f1)
echo "Best scenario is: $BEST_SCENARIO"  # e.g., westG6D1
BEST_RUN_NUM=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f2)
echo "Best run number is: $BEST_RUN_NUM"  # e.g., run48

# Create a folder to perform bootstrapping
mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_boot
# Copy the best run files to the bootstrapping folder
# Using a blerded group name here for easy copying, as the scenario name should be unique across all groups
cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.$BEST_RUN_NUM/$BEST_SCENARIO/$BEST_SCENARIO.pv \
./analysis/DemoInfer/demo_files/${set}_Bestrun_boot/

head -n1 ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.$BEST_RUN_NUM/$BEST_SCENARIO/$BEST_SCENARIO.pv \
    >./analysis/DemoInfer/demo_files/${set}_Bestrun_boot/$BEST_SCENARIO.bootstrap_summary.txt


cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.est \
    ./analysis/DemoInfer/demo_files/${set}_Bestrun_boot/
cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.tpl \
    ./analysis/DemoInfer/demo_files/${set}_Bestrun_boot/
cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.$BEST_RUN_NUM/*_joint*.obs \
    ./analysis/DemoInfer/demo_files/${set}_Bestrun_boot


mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/

grep "^#" ./analysis/DemoInfer/easySFS/${set}.M5F1.syn_inter_LDprune.vcf \
    >./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/$set.vcf.header

grep -v "^#" ./analysis/DemoInfer/easySFS/${set}.M5F1.syn_inter_LDprune.vcf \
    >./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/$set.vcf.body

# Calculate number of sites per block (100 blocks total)
siteCount=$(wc -l < ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/$set.vcf.body)
echo "Total number of sites: $siteCount"

sitesPerBlock=$(($siteCount / 100))
echo "Sites per block: $sitesPerBlock"

mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/blocks
split -l $sitesPerBlock ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/$set.vcf.body \
./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/blocks/$set.sites.

}


# Make bootstrapped SFS files
myfunc_bseasySFS() {
        set="west"
        module load conda
        . activate easySFS


        # Make a new folder for each bootstrapping iteration:
        mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1}

        # Write the header to our new bootstrapped vcf file
        cat ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/$set.vcf.header \
            > ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1}/$set.bs.${1}.vcf

        # Randomly add 100 blocks
        for r in {001..100}; do
            echo "Seed number rep${1}num${r}"

            if [[ -z "$1" ]]; then
                echo "Error: missing repeat number (argument 1)."
                exit 1
            fi

            cat $(shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"seed_rep${1}num${r}" -nosalt </dev/zero 2>/dev/null) -n1 -e ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/blocks/$set.sites.*) \
                >>./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1}/${set}.bs.${1}.vcf
        done

        # Make an SFS from the new bootstrapped file
        python3 ./scripts/easySFS_withoutMultiSFS.py \
            -i ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1}/$set.bs.${1}.vcf \
            -p ./analysis/DemoInfer/easySFS/$set.group.set.txt \
            --order cauc,comD,comP,pyra \
            -a -f \
            --proj 42,28,22,18 \
            -o ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1} \
            --prefix $set
    }




# # Run fastsimcoal2 on bootstrapped SFS files
myfunc_bssbatach() {
        set="west"
        workpath=$(pwd)
        if [[ ! -d "./analysis/DemoInfer" ]]; then
        echo "Error: Script must be run from the project root directory."
        exit 1
        fi
        # Get the best run scenario name and run number
        BEST_SCENARIO=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f1)
        echo "Best scenario is: $BEST_SCENARIO"  # e.g., westG6D1

        BEST_RUN_NUM=$(head -n 2 ./analysis/DemoInfer/demo_files/${set}_summary/Best_run_for_all.summary.txt | tail -n 1 | cut -f2)
        echo "Best run number is: $BEST_RUN_NUM"  # e.g., run48

        re_count=5
        # Parameters
        sim_count=100000
        cpu_per_job=10
        EMC_cycle=40
        # 2026-10-02 rerun: start from best-run .pv (--initValues), -L 40,
        # write to ${set}_Bestrun_block_boot_init; old -L 10 results stay in ${set}_Bestrun_block_boot
        src_dir=./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/${1}
        bs_dir=./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot_init/${1}
        pv_file=./analysis/DemoInfer/demo_files/${set}_Bestrun_boot/$BEST_SCENARIO.pv

        if [[ ! -d $src_dir/fastsimcoal2 ]]; then
            echo "Error: bootstrap SFS not found: $src_dir/fastsimcoal2" >&2; return 1
        fi
        if [[ ! -f $pv_file ]]; then
            echo "Error: best-run .pv not found: $pv_file" >&2; return 1
        fi
        if compgen -G "$bs_dir/$set.run*" > /dev/null; then
            echo "Error: $bs_dir already has run dirs; refusing to overwrite" >&2; return 1
        fi

        # Copy the bootstrap SFS (easySFS output) from the old block_boot dir
        mkdir -p $bs_dir
        cp -rp $src_dir/fastsimcoal2 $bs_dir/
        cp -p $src_dir/datadict.txt $src_dir/easySFS.log $bs_dir/ 2>/dev/null || true

        for i in $(seq -f "%02g" 1 ${re_count}); do
            echo "Prepareing run$i files"
            mkdir -p $bs_dir/$set.run$i
            cp $bs_dir/fastsimcoal2/*_joint*.obs $bs_dir/$set.run$i
            cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.tpl $bs_dir/$set.run$i/$set.tpl
            cp ./analysis/DemoInfer/demo_files/$set/Group_*/$BEST_SCENARIO.est $bs_dir/$set.run$i/$set.est
            cp $pv_file $bs_dir/$set.run$i/$set.pv
            cd $bs_dir/$set.run$i
            cat <<END > bs${1}.run$i.sbatch
#!/bin/bash
#SBATCH -J bs${1}.run$i
#SBATCH -o bs${1}.run$i.%J.out
#SBATCH -e bs${1}.run$i.%J.err
#SBATCH --time=4:00:00

PATH=/shared/projects/pear_snp2/software/fsc28_linux64:$PATH

set -euo pipefail

fsc28 -t ${set}.tpl -e ${set}.est \
    -0 -C 10 -m \
    -n ${sim_count} -L ${EMC_cycle} -M -y 4 \
    --initValues ${set}.pv \
    -q -c ${cpu_per_job} -B 20
END
            sbatch -c ${cpu_per_job} bs${1}.run$i.sbatch
            cd "$workpath"
        done
    }





#########
# Main
#########

## step1
# Test run
# step_fsc_test

## step2
# run fsc
# exit  ## prevent running fsc by accident
# step_fsc_run

## step3
# summary after run fsc, and prepare bootstrap blocks
# step_summary && \
# myfunc_get_boot_bloks

# ## step4
# export -f myfunc_bseasySFS
# bs=($(seq -w 1 50)) # Generate 50 shuffled vcf files
# module load parallel
# parallel -j 25 myfunc_bseasySFS ::: "${bs[@]}"

# ## step5
# export -f myfunc_bssbatach
# bs=($(seq -w 1 50))

# module load parallel
# parallel -j 25 myfunc_bssbatach ::: "${bs[@]}"

# # step6
# mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/summary
# perl ./scripts/fsc-selectbestrun.boot.pl \
#     "./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/*/*.run*/*/*.bestlhoods" \
#     ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot/summary

# module load r/4.3.1
# Rscript ./scripts/boot_intreval.R
