#!/bin/bash
# By Xilong CHEN
# Create date: 2025-06-10
# Contact: chen_xilong@outlook.com

#SBATCH -c 17
#SBATCH --mem=128G
#SBATCH -o twisst5pop.%J.out
#SBATCH -e twisst5pop.%J.err

# source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
# conda activate py310_dtools
# module load parallel/20200822
module load python/3.7
module load parallel/20190322

# chrs=(chr1 chr2 chr3 chr4 chr5 chrX)
chrs=(Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9)

mkdir -p ../output/twisst5pop/twisst5pop.RES

parallel -j 17 \
    "python  ~/work2/software/twisst/twisst.py \
        -t ../output/twisst5pop/genogz2tree/{}.phyml_bionj.5pop.50SNPs.trees.gz \
        -w ../output/twisst5pop/twisst5pop.RES/{}.weights.50SNPs.5pop.csv.gz \
        --outputTopos ../output/twisst5pop/twisst5pop.RES/{}.topologies.50SNPs.5pop.trees \
        -g cauc \
        -g comm_Dessert \
        -g comm_Perry \
        -g pyra \
        -g ussu \
        --outgroup ussu \
        --method complete \
        --groupsFile ../output/twisst5pop/twisst5pop.RES/ussuWest_AB_group.txt" \
    ::: "${chrs[@]}"
