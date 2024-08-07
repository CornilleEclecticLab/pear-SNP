#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Cwd;

# By Xilong CHEN
# Create date: 2022-08-22
# Contact: chen_xilong@outlook.com

# Adapted by Yuqi NIE at 2023-11-29 15:00:03

my $batch= 'test';

my $bin = getcwd();
my $workdir = dirname($bin);
my $input  = "$workdir/input";
my $output = "$workdir/output.test";


my $set_file = "$bin/test_run_Dsuit/s02.set_Dquartets.txt";
my $Dsuite = "$bin/Dsuite/Build/Dsuite";
my $tree_f = "$bin/test_run_Dsuit/s02.tree.txt";
my $chr_f = "$input/s02.chr_list.txt";

my $shdir = "s02.run_Dtest_by_chr.test";
print "Creat folder $shdir\n";
`/bin/rm -rf $shdir`;
`mkdir $shdir`;
`mdkir "$output/$shdir"`

open my $chr_r, "<", $chr_f or die "Could not open $chr_f:$!\n";
my @chrs;
while (my $line = <$chr_r>) {
    chomp $line;
    push @chrs, $line;
}
close $chr_f;


for my $chr (@chrs) {
    my $input_vcfgz =
      "$bin/test_run_Dsuit/$chr.choosing.test.vcf.gz";
    my $out_prefix = "$output/$shdir/$chr";
    my $job = "$shdir.$chr";
    my $sh  = << "SH_END";
#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err

source $bin/s00.load_environment.sh

$Dsuite Dtrios \\
    --KS-test-for-homoplasy \\
    -o $out_prefix \\
    -t $tree_f \\
    $input_vcfgz \\
    $set_file

SH_END

    open my $ou, ">", "$shdir/$job.sh";

    print $ou $sh;

}
