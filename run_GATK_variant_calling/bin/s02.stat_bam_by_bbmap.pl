#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/14/2020 16:05
# Contact: chen_xilong@outlook.com

my $shdir = "s02.stat_bam_by_bbmap";
`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;

my $bbm     = "pileup.sh";
my $workdir = "~/work/pear/run_GATK_variant_calling";
my $output  = "$workdir/output/$shdir";

open my $fh, "<", "bam_list.txt";

my @bamlist;

# while (<$fh>) {
#     chomp;
#     push @bamlist, $_;
# }

# my @chunks;
# my $size = 20;
# push @chunks, [ splice( @bamlist, 0, $size ) ] while @bamlist;
# my $index = 1;
# foreach my $ref (@chunks) {
#     open my $ou, ">", "$shdir/$index.sh";
#     for (@$ref) {
#         my @aa = fileparse $_;
#         my $d  = basename $aa[1];
#         print $ou
#             "$bbm -Xmx2g in=$_ out=$output/$d.stats > $output/$d.bamreads 2>&1 && \\\n";
#     }
#     $index += 1;
# }

while (<$fh>) {
    chomp;
    my @aa  = fileparse $_;
    my $d   = basename $aa[1];
    my $job = "$d.bbmap";
    my $sh  = << "SH_END";
#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err
module purge

module load bioinfo/samtools-1.9
module load bioinfo/bbmap_38.95

$bbm -Xmx2g in=$_ out=$output/$d.stats > $output/$d.bamreads 2>&1

SH_END

    open my $ou, ">", "$shdir/$job.sh";
    print $ou $sh;
}
