#!/bin/bash

#SBATCH -J s01.combine
#SBATCH -e s01.b1.%j.err
#SBATCH -o s01.b1.%j.out

# 2023/01/25 modify this for sbatch running
# 

module purge 
module load bcftools/1.14

list=$1 
prefix=$2
combine_vcf=$list
test_vcf="$prefix.test.vcf.gz"
sum_number_txt="$prefix.sum_number.txt"
ind_miss_txt="$prefix.ind_miss.txt"
test_rand_vcf="$prefix.test_rand.vcf"
test_rand_vcfgz="$prefix.test_rand.vcf.gz"
GotTest="False"


bcftools concat -f $list -O z4 -o $test_rand_vcfgz

tabix -p vcf $test_rand_vcfgz


# count genotype number 
bcftools view -H $test_rand_vcfgz \
| wc -l \
> $sum_number_txt


# count SNP number for each individual
bcftools stats -s - $test_rand_vcfgz \
| grep -E ^PSC \
| cut -f3,14 \
> $ind_miss_txt

