#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-03-22
# Contact: chen_xilong@outlook.com

my $shdir = "sx2.concat.count";
`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;
print "Creat folder $shdir\n";

my $workpath =
  "/work/xchen/XCHEN_apple_genome/analysis_2022_Mar/run_GATK_variant_calling";
my $output = "$workpath/output";
my $input  = "$workpath/input";

print "Creat output folder $output/$shdir\n";
`mkdir $output/$shdir`;

my @chrs;
for my $i ( 1 .. 17 ) {
    my $j = sprintf( "%#02s", $i );
    push @chrs, "Chr$j";
}

for my $chr (@chrs) {
    my $vcf    = "/work/xchen/XCHEN_apple_genome/analysis_2022_Mar/input/concat_vcf/$chr.filtered_pixy_concat.vcf.gz";
    my $outcount = "$output/$shdir/$chr.filtered_pixy_concat.count.txt";
    my $job    = "$chr.count";
    my $sh     = <<"SH_END";
#!/bin/bash

#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err

module purge
module load bioinfo/bcftools-1.9
module load bioinfo/tabix-0.2.5

bcftools view -H $vcf | wc -l > $outcount

SH_END
    open my $ou, ">", "$shdir/$job.sh";
    print $ou $sh;
}

