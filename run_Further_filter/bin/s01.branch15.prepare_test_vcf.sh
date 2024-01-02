#!/bin/bash

#SBATCH -J s01.prepare
#SBATCH -e s01.slurm.%J.err
#SBATCH -o s01.slurm.%J.out
#SBATCH -c 8

# 2023/01/25 modify this for sbatch running
# 2023-07-30 modify this for bcftools concat
# 2023-11-30 15:18:24 This scipt is used to get random chosed data from vcf, which could be used as three mode:
    # 1. If you have vcf for each scaffold or chromosome, you want to combine them, and then randome chose data,
        # run this scirpt using paramater $1 and $2.
    
    # 2. If you have vcf for each scaffold or chromosome, but you want to randome chose data from each vcf, and then combine them,
        # run this script using parameter $1, $2 and $3, $1 coubld be any words, as it will be ignored
    
    # 2. If you only have one combined vcf.gz for whole genome, 
        # run this script with $1, $2, and $3, $1 could be any words, as it will be ignored
        # and $3 is the path of the combined vcf.gz file

module purge 
module load bcftools/1.14
 


usage="bash $0 <vcf_list_file> <output_prefix> [one_vcf_file_path]\n"

if [ $# -lt 2 ]
then
    echo -e "Usage:\n"
    echo -e $usage
    echo -e "When one vcf file path is given as the third paramater, the first parameter will be ignored\n"
    exit 1
fi


list=$1          # 2023-11-30 16:04:16 The $1 parameter is necessary
prefix=$2        # 2023-11-30 16:04:16 
combine_vcf=$3   # 2023-11-30 16:04:16 However, if the $3 parameters is not empty, the $1 parameter will be ignored
test_vcf=$4      # 2023-11-30 16:09:02 If the $4 parameters is not empty, the first and the third parameters will be ignored

GotTestFile="False"
GotOneCombinedFile="False"

if [ "$combine_vcf" != "" ]; then
    GotOneCombinedFile="True"
fi

if [ "$test_vcf" != "" ]; then
    GotTestFile="True"
fi


test_vcf="$prefix.test.vcf.gz"
sum_number_txt="$prefix.sum_number.txt"
ind_miss_txt="$prefix.ind_miss.txt"
test_rand_vcf="$prefix.test_rand.vcf"
test_rand_vcfgz="$prefix.test_rand.vcf.gz"


# The default parameter '-O z' using compress level 9 that is very slow


### According to the test on 2023-07-31, the following is convient but about 90% slower than seperately run filter on each chromosome vcf.gz file
### Updated this part on 2023-1204
if [ "$GotTestFile" == "False" ] && [ "$GotOneCombinedFile" == "False" ]; then
    bcftools concat --threads 8 -f $list \
    | bcftools view -m2 -M2 -v snps -O z4 -o $test_vcf --threads 8 
    # | bcftools filter -e 'F_MISSING > 0.2' -O z4 -o $test_vcf --threads 8
    # | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 -o $test_vcf --threads 8
    GotTestFile="True"
    echo "GotTestFile is $GotTestFile"
fi


### Alternatively, run filter on each chromosome vcf.gz file seperately, without combine them, or you have already combined them
if [ "$GotTestFile" == "False" ] && [ "$GotOneCombinedFile" == "True" ]; then
    bcftools view $combine_vcf -m2 -M2 -v snps -O z4 -o $test_vcf --threads 8 
    # |bcftools filter -e 'F_MISSING > 0.2' -O z4 -o $test_vcf --threads 8
    #|bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 -o $test_vcf --threads 8
    GotTestFile="True"
fi


tabix -p vcf $test_vcf


# choose random 1% site for DP plot and clone test and missing test
perl ./SNP_rand_choose_gz.pl 100 $test_vcf $test_rand_vcf


bgzip $test_rand_vcf
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


# remove the combine vcf
#rm $combine_vcf
#rm $test_vcf
