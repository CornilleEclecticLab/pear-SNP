#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 4/11/2020 19:05
# Contact: chen_xilong@outlook.com

# Edited by Yuqi NIE
# For IFB. And Copy meanQ files to a folder for CLUMPAK
# Edit Date: 03/08/2023 10:45

my $FastStructure_py = "/shared/ifbstor1/software/miniconda/envs/faststructure-1.0/bin/structure.py";
my $usage
    = "Usage: perl $0 -b <bedfile name(not include .bed)> -minK <minK> -maxK <maxK> -r <repeat number> -o <output_dir>
    \nNote1: You may need to modify the directory path of structure.py in the script (now the path is $FastStructure_py).\n\n";



# print $usage;

my $bedfile;
my $output_dir = ".";
my $minK       = 2;
my $maxK       = 19;
my $repeat     = 2;

## parse options from @ARGV
GetOptions(
    'b=s'    => \$bedfile,
    'o=s'    => \$output_dir,
    'minK=i' => \$minK,
    'maxK=i' => \$maxK,
    'r=i'    => \$repeat
) or die $usage;

die("ERROR: bedfile path must be specified.\n$usage") unless defined $bedfile;

for my $k ( $minK .. $maxK ) {
    `mkdir $output_dir/K$k`;
    for ( 1 .. $repeat ) {
        my $seed  = int(rand(100000));
        my $job = "K$k.run$_.sh";
        open my $ou, ">", "$output_dir/K$k/run$_.sh";
            my $fssh = <<"END_FSSH";
#!/bin/bash

#SBATCH -J $job
#SBATCH -o $job.%J.out
#SBATCH -e $job.%J.err

module purge
module load faststructure

python $FastStructure_py -K $k --cv 10 --full --seed=$seed --input=$bedfile --output=run$_

END_FSSH
        # module load bioinfo/fastStructure-1.0
        #
        print $ou "$fssh\n";
        close $ou;
    }
}


open my $ou1, ">", "00launch_submits.sh" or die;
print $ou1 "
cd $output_dir
for rep in `seq $minK $maxK`
do
        cd ./K\$rep
        ls ./ | xargs -I {} sbatch -c 1 --mem=1G {}
        cd ..
done";

open my $ou2, ">", "01.zip_for_clumpak.sh" or die;
print $ou2 "
cd $output_dir
mkdir result
cp ./K*/*.meanQ ./result

cd result
for i in `seq $minK $maxK`
do
    zip -q K\${i}.zip *.\${i}.meanQ
done

zip -q ../input_clumpak.zip K*.zip

rm -rf K*.zip
cd $output_dir
rm -rf ./result
"
