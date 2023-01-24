#!/bin/bash


usage="bash $0 <vcf_list_file> <output_prefix>\n"

if [ "$#" -ne 2 ]
then
    echo -e $usage
    exit 1
fi

list="$1";
combine_vcf="$2.combine.vcf.gz"
test_vcf="$2.test.vcf.gz"
sum_number_txt="$2.sum_number.txt"
ind_miss_txt="$2.ind_miss.txt"
test_rand_vcf="$2.test_rand.vcf"
test_rand_vcfgz="$2.test_rand.vcf.gz"

# if you have vcf for each scaffold or chromosome, combine them
# vcf.file.list contain the path of each vcf.gz file
bcftools concat -f $list -O z \
 -o $combine_vcf

# if you only have one vcf.gz for whole genome, you can start here
# ectract all the SNP and filter genotype missing 20% 
bcftools view -v snps \
$combine_vcf \
| bcftools filter -e 'F_MISSING > 0.2' -O z \
-o $test_vcf

tabix -p vcf $test_vcf

# choose random 1% site for DP plot and clone test and missing test
perl ./SNP_rand_choose_gz.pl 100 $test_vcf $test_rand_vcf

bgzip $test_rand_vcf
tabix -p vcf $test_rand_vcfgz

# count genotype number 
bcftools view -H \
$test_rand_vcfgz \
| wc -l > $sum_number_txt

# count SNP number for each individual
bcftools stats -s - $test_rand_vcfgz \
| grep -E ^PSC \
| cut -f3,14 \
> $ind_miss_txt


# remove the combine vcf
#rm $combine_vcf
#rm $test_vcf
