#!/bin/bash
# By Xilong CHEN
# Create date: 2025-08-17
# Contact: chen_xilong@outlook.com

# 2025-08-21 updated for comm ref

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
    mkdir ./analysis/DemoInfer/PrepareVcf

    bcftools filter -e 'F_MISSING > 0.2' \
        ./data/fsc_input_250816/pear_Jul2024_ref_comm.Europe_noAdmix.Combine_chr.maf005.vcf.gz -Oz \
        -o ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz
    tabix -p vcf ./analysis/DemoInfer/PrepareVcf/west.M5F1.gz

    # VEP
    mkdir ./analysis/DemoInfer/PrepareVcf/syn_inter
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

    mkdir ./analysis/DemoInfer/PrepareVcf/plink_LDprune
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
    mkdir ./analysis/DemoInfer/easySFS
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
    mkdir ./analysis/DemoInfer/demo_files
    mkdir ./analysis/DemoInfer/demo_files/west
    # a test run for west folder, the script will only run once with
    # 100 simulation, and make the plot. It's better test them one by one,
    # The running time is very short, so it won't take up time, but the
    # output/error and the demo plots can help identify errors in tpl and est.
    # I modified the ParFileViewer.r script to ParFileViewer.logname.r here
    # to scale the y-axis in log10 and accept the pop names parameter.
    # If your model fits better with the original version, please use the original script.

    set="west"
    workpath="/home/xchen/work/PearPop_Project_Yuqi"
    Group="Group_G0"
    prefix_list=("westG0D1" "westG0D2" "westG0D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG1D1" "westG1D2" "westG1D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG2D1" "westG2D2" "westG2D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG3D1" "westG3D2" "westG3D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG4D1" "westG4D2" "westG4D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG5D1" "westG5D2" "westG5D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG6D1" "westG6D2" "westG6D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG7D1" "westG7D2" "westG7D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG8D1" "westG8D2" "westG8D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    prefix_list=("westG9D1" "westG9D2" "westG9D3")
    for prefix in "${prefix_list[@]}"; do {
        mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.test_run
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
    # workpath="/home/xchen/work/PearPop_Project_Yuqi"
    workpath="/scratch/yn2515/pear_snp_tip/run_fastsimcoal2"
    # Group_list=("Group_G0" "Group_G1" "Group_G2" "Group_G3")

    G_list=("G0" "G1" "G2" "G3" "G4" "G5" "G6" "G7" "G8" "G9")
    for g in "${G_list[@]}"; do
        prefix_list=()
        Group="Group_${g}"
        for d in 1 2 3; do
            prefix_list+=("west${g}D${d}")
        done
        echo $Group
        echo "${prefix_list[@]}"

        for prefix in "${prefix_list[@]}"; do
            for i in $(seq -w 1 50); do
                echo "Prepareing run$i files"
                mkdir ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i
                cp ./analysis/DemoInfer/easySFS/fastsimcoal2/*_joint*.obs \
                    ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i
                cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.tpl \
                    ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i
                cp ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.est \
                    ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i

                cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i
                for file in *.obs; do
                    newname=$(echo "$file" | sed "s/$set/$prefix/g")
                    mv "$file" "$newname"
                done
                cd $workpath || exit
            done
        done
        # submit jobs
        for prefix in "${prefix_list[@]}"; do
            for i in $(seq -w 1 50); do
                cd ./analysis/DemoInfer/demo_files/west/$Group/${prefix}.run$i
                cat <<END >$prefix.run$i.sbatch
#!/bin/bash
#SBATCH -J $prefix.run$i
#SBATCH -o $prefix.run$i.%J.out
#SBATCH -e $prefix.run$i.%J.err

PATH=/home/yn2515/software/fsc28_linux64:$PATH

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
    mkdir ./analysis/DemoInfer/demo_files/west_summary
    perl ./scripts/fsc-selectbestrun.pl \
        "./analysis/DemoInfer/demo_files/west/*/*.run*/*/*.bestlhoods" \
        ./analysis/DemoInfer/demo_files/west_summary
    cat ./analysis/DemoInfer/demo_files/west_summary/Best_run_for_all.summary.txt
    # Scenario        run_nu  MaxEstLhood     AIC     deltaL  best_of_what
    # westG6D1        run48   -1833886.714    8445380.41979247        96438.9639999999        best_of_MaxEstLhood
    # westG6D1        run48   -1833886.714    8445380.41979247        96438.9639999999        best_of_AIC
    mkdir ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file

    cp ./analysis/DemoInfer/demo_files/west/Group_G6/westG6D1.run48/westG6D1/*.bestlhoods \
        ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file/
    cp ./analysis/DemoInfer/demo_files/west/Group_G6/westG6D1.run48/westG6D1/*.par \
        ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file/

    workpath="/home/xchen/work/PearPop_Project_Yuqi"
    cd ./analysis/DemoInfer/demo_files/west_Bestrun_bestlhoods_file/
    Rscript $workpath/scripts/ParFileViewer.logname.r \
        westG6D1_maxL.par \
        cauc,comD,comP,pyra
    cd $workpath

# compare the scenarios
    Rscript ./scripts/fsc_MaxEsthood_boxplot.R
    # Rscript ./scripts/fsc_best_confidence_interval.R

}

# step_fsc_run
