#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 3/9/2021 23:27
# Contact: chen_xilong@outlook.com

# before this step you have to know DP cutoff and remove.list of high missing 
# indvidual and clone individual

my $shdir = "s05.filter_concat";

`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;

my @Chrs = qw(
GWHBAOS00000076
GWHBAOS00000158
GWHBAOS00000386
GWHBAOS00000381
GWHBAOS00000224
GWHBAOS00000085
GWHBAOS00000163
GWHBAOS00000172
GWHBAOS00000128
GWHBAOS00000425
GWHBAOS00000352
GWHBAOS00000424
GWHBAOS00000274
GWHBAOS00000099
GWHBAOS00000356
GWHBAOS00000365
GWHBAOS00000335
scaffolds.list
);

my $workdir = "/home/zruilin/work/Peach/run_Further_filter";
my $vcfdir = "$workdir/input/vcfs";
my $outdir = "$workdir/output/$shdir";

for my $chr (@Chrs) {
    my $INPUT = "$vcfdir/$chr.rmQFI.combine.vcf.gz";
    my $OU1   = "$outdir/$chr.filtered_invariant.vcf.gz";
    my $OU2   = "$outdir/$chr.filtered_variant.tmp.vcf.gz";
    my $OU3   = "$outdir/$chr.filtered_concat.vcf.gz";
    my $variant_vcf = "$outdir/$chr.filtered_variant.vcf.gz";
    my $outcount_concat = "$outdir/$chr.filtered_concat.count.txt";
    my $outcount_variant = "$outdir/$chr.filtered_variant.count.txt";
    my $job = "$chr.s05";
    my $p     = <<"S01END";
#!/bin/bash

#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err

module load bioinfo/bcftools-1.9
module load bioinfo/tabix-0.2.5

# extract the invariants and do filtering
bcftools filter -S . -e 'FMT/DP<5 | FMT/RGQ<20 | FMT/DP>100' \\
$INPUT \\
| bcftools filter -i 'ALT="."' \\
| bcftools filter -e 'F_MISSING > 0.2' -O z \\
-o $OU1

# extract the variants and do filtering
bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20 | FMT/DP>100' \\
$INPUT \\
| bcftools filter --SnpGap 10 \\
| bcftools view -m2 -M2 -v snps \\
| bcftools filter -e 'F_MISSING > 0.2' -O z \\
-o $OU2

tabix -p vcf $OU1
tabix -p vcf $OU2

# creat a concat vcf for all the variants and invariants
bcftools concat \\
--allow-overlaps \\
$OU1 $OU2 \\
-O z -o $OU3

tabix -p vcf $OU3

# count the site number of concat vcf files
bcftools view -H $OU3 | wc -l > $outcount_concat


# as we have the concat vcf, we can further filter the AN=AC and AC=0 sites in variants
# also count the SNP number

bcftools filter -e 'AC==0 || AC==AN' -O z -o $variant_vcf $OU2 && \\
tabix -p vcf $variant_vcf  && \\
bcftools view -H $variant_vcf | wc -l > $outcount_variant

# remove two temp vcf files
rm $OU1
rm $OU2

S01END
    open my $ou, ">", "./$shdir/$job.sh";
    print $ou $p;
}
