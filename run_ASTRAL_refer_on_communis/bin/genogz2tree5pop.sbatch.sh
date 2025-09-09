#!/bin/bash
# By Xilong CHEN
# Create date: 2025-06-10
# Contact: chen_xilong@outlook.com
# Modified by Yuqi on 2025-09-01


#SBATCH --job-name=genogz2tree5pop
#SBATCH -e genogz2tree5pop.%j.err
#SBATCH -o genogz2tree5pop.%j.out
#SBATCH -c 68
#SBATCH --mem=512G

# source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
# conda activate py310_dtools
# module load parallel/20200822
module load parallel/20190322
module load python/3.7
module load phyml/3.3.20190909


# chrs=(Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9)
readarray -t chrs < ../input/s01.chr_list.txt

mkdir -p ../output/twisst5pop/genogz2tree


export PYTHONPATH=./genomics_general.NOGIT.:$PYTHONPATH

parallel -j 17 \
    "python ./genomics_general.NOGIT./phylo/phyml_sliding_windows.py \
    -T 4 \
    -g ../output/vcfgz2genogz/{}.ussu.west.noAdmix.phased.geno.gz \
    --prefix ../output/twisst5pop/genogz2tree/{}.phyml_bionj.5pop.50SNPs \
    -w 50 \
    --windType sites \
    --model GTR \
    --outgroup SRR7135510,WPU1,WPU10,WPU11,WPU12,WPU14,WPU15,WPU3,WPU4,WPU5,WPU7 \
    --optimise n" \
::: "${chrs[@]}"
