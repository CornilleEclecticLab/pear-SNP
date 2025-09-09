#!/bin/bash
# By Xilong CHEN
# Create date: 2025-06-10
# Contact: chen_xilong@outlook.com
# Modified by Yuqi on 2025-09-01


#SBATCH -c 17
#SBATCH --mem=128G
#SBATCH -e vcfgz2genogz5pop.%j.err
#SBATCH -o vcfgz2genogz5pop.%j.out

module purge
# module load parallel/20200822
#source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
# conda activate py310_dtools

# module load parallel/20200822
module load parallel/20190322
module load python/3.7


# chrs=(chr1 chr2 chr3 chr4 chr5 chrX)
chrs=(Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9)

mkdir -p ../output/vcfgz2genogz

export PYTHONPATH=./genomics_general.NOGIT.:$PYTHONPATH

parallel -j 17 \
'sh -c "python ./genomics_general.NOGIT./VCF_processing/parseVCF.py \
    -i ../output/vcfgz/{}.west.noAdmix.phased.vcf.gz |
    gzip > ../output/vcfgz2genogz/{}.ussu.west.noAdmix.phased.geno.gz"' \
::: "${chrs[@]}"


