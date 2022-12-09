#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-03-31
# Contact: chen_xilong@outlook.com
# Adapted by Nie Yuqi on 2022-10-19

my $workdir = "/work/zruilin/Peach/QC";
my $output  = "$workdir/output";
my $input   = "$workdir/input";

my $shdir = "s01.raw_qc";
print "creat $shdir
use 2 cpu\n";
`rm -rf $shdir`;
`mkdir $shdir`;

open my $in, "<", "$input/raw_data.list" or die;
my %Dnu;
while (<$in>) {
    chomp;
    my @aa  = fileparse $_;
    my $dnu = basename $aa[1];
    print "$dnu\n";
    push @{ $Dnu{$dnu} }, $_;
}

 print Dumper \%Dnu;

for my $d ( sort keys %Dnu ) {
    print "$d\n";
    my $job = "$d.$shdir";
    open my $ou, ">", "$shdir/$job.sh";
    my @f = sort @{ $Dnu{$d} };
    print $ou "#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err
module purge
  
module load bioinfo/FastQC_v0.11.7
mkdir $output/$shdir/$d
";
    for my $in_f (@f) {
        print $ou "fastqc -t 2 -f fastq -o $output/$shdir/$d $in_f &&\\ \n";
    }
    print $ou "echo done\n";
}
