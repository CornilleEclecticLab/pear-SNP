#!/bin/bash
# By Xilong CHEN
# Create date: 2025-06-09
# Contact: chen_xilong@outlook.com
# Modified by Yuqi on 2025-09-01
#SBATCH --time=24:00:00
#SBATCH -c 18
#SBATCH -e bcf2vcfgz5pop.sbatch.%j.err
#SBATCH -o bcf2vcfgz5pop.sbatch.%j.out



module purge
# module load parallel/20200822
module load parallel/20190322

#source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
# conda activate Aphid
module load bcftools/1.14

# chrs=(chr1 chr2 chr3 chr4 chr5 chrX)
chrs=(Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9)

mkdir -p ../output/vcfgz

parallel -j 17 \
"bcftools view -S ../input/ussu.west.list.txt \
    ../input/{}.phased.reheader.vcf.gz \
    -Oz \
    -o ../output/vcfgz/{}.west.noAdmix.phased.vcf.gz &&
    tabix -p vcf ../output/vcfgz/{}.west.noAdmix.phased.vcf.gz" \
::: "${chrs[@]}"

