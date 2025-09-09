#!/bin/bash

jobID=1

step_twisst5pop() {
    # 4 pure west pop, with ussu as 5th pop and outgroup
    mkdir -p ../output/twisst5pop

    mkdir -p ../output/twisst5pop/vcfgz

    jobID=$(sbatch --parsable  -J bcf2vcfgz5pop -c 6 -o ./sbatch_log/bcf2vcfgz5pop.log \
        ./bcf2vcfgz5pop.sbatch.sh)

    mkdir ../output/twisst5pop/vcfgz2genogz
    jobID=$(sbatch --parsable --dependency=afterok:$jobID -J vcfgz2genogz5pop -c 17 -o ./sbatch_log/vcfgz2genogz5pop.log \
        ./vcfgz2genogz5pop.sbatch.sh)

    mkdir ../output/twisst5pop/genogz2tree
    jobID=$(sbatch --parsable --dependency=afterok:$jobID -J genogz2tree5pop -c 24 -o ./sbatch_log/genogz2tree5pop.log \
        ./genogz2tree5pop.sbatch.sh )


    mkdir -p ../output/twisst5pop/twisst5pop.RES
    cat ../input/s01.west.sample_tab_population.txt |
        awk '{print $1"_A""\t"$2}' \
            > ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt
    cat ../input/s01.west.sample_tab_population.txt |
        awk '{print $1"_B""\t"$2}' \
            >> ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt
    awk '{print $1"_A""\tussu"}' ../input/ussu.list.txt \
        >> ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt
    awk '{print $1"_B""\tussu"}' ../input/ussu.list.txt \
        >> ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt

    jobID=$(sbatch --parsable --dependency=afterok:$jobID -J twisst5pop -c 6 -o ./sbatch_log/twisst5pop.log \
        ./twisst5pop.sbatch.sh)

}

step_ASTER5pop() {
    mkdir -p ../output/ASTER.5pop

    mkdir -p ../output/ASTER.5pop/filter_tree

    module load parallel
    module load bedtools

    #chrs=(Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9)
    readarray -t chrs < ../input/s01.chr_list.txt

    parallel -j 17 \
        "perl ./rm_masked_tree.pl \
        ../input/Pyrus_communis.genmap.mask.bed \
        ../output/twisst5pop/genogz2tree/{}.phyml_bionj.5pop.50SNPs.data.tsv \
        ../output/twisst5pop/genogz2tree/{}.phyml_bionj.5pop.50SNPs.trees.gz \
        0.2 \
    > ../output/ASTER.5pop/filter_tree/{}.phyml_bionj.5pop.50SNPs.trees" \
        ::: "${chrs[@]}"


    # prepare group
    mkdir -p ../output/ASTER.5pop/astral4.5pop.RES
    cp ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt \
       ../output/ASTER.5pop/astral4.5pop.RES/
    # prepare trees
    cat ../output/ASTER.5pop/filter_tree/Chr*.phyml_bionj.5pop.50SNPs.trees \
        > ../output/ASTER.5pop/astral4.5pop.RES/auto.5pop.trees

    sbatch -J astral4.auto.5pop \
        ./astral4.auto.5pop.sbatch.sh
}

step_twisst5pop

Wait for the previous job to complete before running the next step
while squeue -j $jobID 2>/dev/null | grep -q $jobID; do
    sleep 10
done

step_ASTER5pop
