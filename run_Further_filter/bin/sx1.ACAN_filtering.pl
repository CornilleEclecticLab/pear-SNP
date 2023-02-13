#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/20/2020 18:03
# Contact: chen_xilong@outlook.com

my $shdir = "sx1.ACAN_filtering";
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
    my $vcf    = "$input/variant_vcf/$chr.filtered_pixy_variant.vcf.gz";
    my $outvcf = "$output/$shdir/$chr.filtered_variant.vcf.gz";
    my $outcount = "$output/$shdir/$chr.filtered_variant.count.txt";
    my $job    = "$chr.ACAN";
    my $sh     = <<"SH_END";
#!/bin/bash

#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err

module purge
module load bioinfo/bcftools-1.9
module load bioinfo/tabix-0.2.5

bcftools filter -e 'AC==0 || AC==AN' -O z -o $outvcf $vcf && \\
tabix -p vcf $outvcf  && \\
bcftools view -H $outvcf | wc -l > $outcount

SH_END
    open my $ou, ">", "$shdir/$job.sh";
    print $ou $sh;
}
