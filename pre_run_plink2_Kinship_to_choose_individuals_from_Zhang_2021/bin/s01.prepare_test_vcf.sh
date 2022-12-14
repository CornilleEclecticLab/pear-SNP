#!/usr/bin/env bash


# @File    :   s01.prepare_test_vcf.sh
# @Version :   1.0
# @Author  :   Nie Yuqi
# @Email   :   nieyuqi.cn@gmail.com
# @Time    :   2022/12/06 14:55:09
#Description:
    # 
module purge
module load bioinfo/bcftools-1.16
module load bioinfo/tabix-0.2.5

# Record software version
info='info.software_version.txt'
bcftools version >> $info
echo -e "\n tabix Version: 0.2.5(r964)\n" >> $info
perl -version >> $info

usage="Usage: \n  bash $0 <vcf_file> <output_prefix> \n\nExample: \n  bash $0 ../input/genotype_GATK.vcf.gz ../output/Pear.Zhang_2021 \n"


if [ "$#" -ne 2 ]
then 
    echo -e $usage
    exit 1
fi


vcf="$1"
test_vcf="$2.test.vcf.gz"
sum_before_number_txt="$1.sum_number.txt"
ind_before_miss_txt="$1.ind_miss.txt"
sum_number_txt="$2.sum_number.txt"
ind_miss_txt="$2.ind_miss.txt"
test_rand_vcf="$2.test_rand.vcf"
test_rand_vcfgz="$2.test_rand.vcf.gz"


# ectract all the SNP and filter genotype missing 20%
bcftools view -v snps \
$vcf \
| bcftools filter -e 'F_MISSING > 0.2' -O z \
-o $test_vcf


tabix -p vcf $test_vcf


# choose random 1% site for DP plot and clone test and missing test
perl ./bin.SNP_rand_choose_gz.pl 100 $test_vcf $test_rand_vcf 


bgzip -f $test_rand_vcf 
tabix -p vcf $test_rand_vcfgz


# count genotype number 
bcftools view -H \
$test_rand_vcfgz \
| wc -l > $sum_number_txt

# bcftools view -H \
# $vcf \
# | wc -l > $sum_before_number_txt


# count SNP number for each individual
bcftools stats -s - $test_rand_vcfgz \
| grep -E ^PSC \
| cut -f3,14 \
> $ind_miss_txt

# bcftools stats -s - $vcf \
# | grep -E ^PSC \
# | cut -f3,14 \
# > $ind_before_miss_txt


# remove the combine vcf
## rm $test_vcf