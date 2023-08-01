#!/bin/bash

#SBATCH -J s01.prepare
#SBATCH -e s01.%J.err
#SBATCH -o s01.%J.out

# 2023/01/25 modify this for sbatch running
# 2023-07-30 modify this for bcftools concat

module purge 
module load bcftools/1.14
 


# usage="bash $0 <vcf_list_file> <output_prefix>\n"

# if [ "$#" -ne 2 ]
# then
#    echo -e $usage
#    exit 1
# fi



prefix='pear_Jul2023' #$2
list=$4
combine_vcf=$3
prefix=$2
test_vcf=$1
GotOneCombinedFile="False"
GotTestFile="False"

if [ "$combine_vcf" != "" ]; then
    GotOneCombinedFile="True"
fi

if [ "$test_vcf" != "" ]; then
    GotTestFile="True"
fi

# list='/shared/ifbstor1/projects/pear_snp3/pear/run_Further_filter/input/s01.branch8.input.pear_Jul2023.filter_passed_sites_vcf.file_list.txt'
# list='/shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/input/s01.input.whole_pear_filter_passed_sites_vcf.file_list.txt'
# prefix='whole_pear'

test_vcf="$prefix.test.vcf.gz"
# test_vcf='/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/branch8.pear_Jul2023/pear_Jul2023.Combine_Chr.geno20_maf005.vcf.gz' #2023-07-29Changed from test_vcf="$prefix.test.vcf.gz"
sum_number_txt="$prefix.sum_number.txt"
ind_miss_txt="$prefix.ind_miss.txt"
test_rand_vcf="$prefix.test_rand.vcf"
test_rand_vcfgz="$prefix.test_rand.vcf.gz"




# If you have vcf for each scaffold or chromosome, combine them
# vcf.file.list contain the path of each vcf.gz file
# The default parameter '-O z' using compress level 9 that is very slow


### I have done the fllowing steps in /shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/branch8.pear_Jul2023/pear_Jul2023.Combine_Chr.geno20_maf005.vcf.gz' # test_vcf="$prefix.test.vcf.gz"

# if [ "$GotTest" == "False" ]; then
    # bcftools concat -f $list \
    # | bcftools view -v snps \
    # | bcftools filter -e 'F_MISSING > 0.2' -O z4 -o $test_vcf 
    # $GotTest="True"
# fi

### According to the test on 2023-07-31, the following is convient but about 90% slower than seperately run filter on each chromosome vcf.gz file
### Updated this part on 2023-07-30
if [ "$GotTestFile" == "False" ] && [ "$GotOneCombinedFile" == "False" ]; then
    bcftools concat --threads 8 -f $list \
    | bcftools view -m2 -M2 -v snps --threads 8 \
    | bcftools filter -e 'F_MISSING > 0.2' -O z4 -o $test_vcf --threads 8
    # | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 -o $test_vcf --threads 8
    GotTestFile="True"
    echo "GotTestFile is $GotTestFile"
fi


# If you only have one vcf.gz for whole genome, you can change the last
# paragraph code as annotation and start here
# ectract all the SNP and filter genotype missing 20% 

if [ "$GotTestFile" == "False" ] && [ "$GotOneCombinedFile" == "True" ]; then
    bcftools view -m2 -M2 -v snps --threads 8 $combine_vcf \
    |bcftools filter -e 'F_MISSING > 0.2' -O z4 -o $test_vcf --threads 8
    #|bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05 || AC==0 || AC==AN' -O z4 -o $test_vcf --threads 8
    GotTestFile="True"
fi


tabix -p vcf $test_vcf


# choose random 0.1% site for DP plot and clone test and missing test
perl ./SNP_rand_choose_gz.pl 1000 $test_vcf $test_rand_vcf


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
