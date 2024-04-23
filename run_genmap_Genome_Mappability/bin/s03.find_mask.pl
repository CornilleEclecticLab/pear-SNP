#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw(sum);

# By Xilong CHEN
# Create date: 2021-12-24
# Contact: chen_xilong@outlook.com

# This had two ISSUEs:
#  1. The start position of bed file should be 0-based.
#  2. This script calculate the scores from 0 and 0+step seperatelly. 
# When the scaffold smaller than one window size, the sencond parts (0+steps to end) no needs calculated. 

# Update: 2024-04-23 17:23:24
# Fix the issue 1.

my $infile = shift @ARGV;
my $seed = int(rand() * 1000000);
my $tmp0 = "$seed.NOGIT.tmp0";
my $tmp1 = "$seed.NOGIT.tmp1";
my $output = shift @ARGV;

open my $in,  '<', "$infile" or die;
open my $ou0, ">", "$tmp0"   or die;
open my $ou1, ">", "$tmp1"   or die;
my $window = 100000;
my $step = 50000;
$/ = ">";
<$in>;

while (<$in>) {
    chomp;
    my ( $chr, $seq ) = split /\n/, $_, 2;
    my @values = split / /, $seq;
    for ( my $i = 0 ; $i < ( $@values ) ; $i += $window ) {
        my $j     = $i + $window;
        $j = $@values if $j > $@values;
        my @slice = @values[ $i .. $j - 1];
        my $avg   = sum(@slice) / @slice;
        my $start = $i;
        my $end   = $j;
        if ( $avg < 0.9 ) {
            print $ou0 "$chr\t$start\t$end\t$avg\n";
        }
    }
        for ( my $i = $step ; $i < ( $@values ) ; $i += $window ) {
        my $j     = $i + $window;
        $j = $@values if $j > $@values;
        my @slice = @values[ $i .. $j - 1 ];
        my $avg   = sum(@slice) / @slice;
        my $start = $i;
        my $end   = $j;
        if ( $avg < 0.9 ) {
            print $ou1 "$chr\t$start\t$end\t$avg\n";
        }
    }
}

`cat $tmp0 $tmp1 > $output.bed`;
`sort -k1,1 -k2,2n $output.bed > $output.100k.sorted.bed`;
`module load bedtools && \
 bedtools merge -i $output.100k.sorted.bed > $output.100k.sorted.merge.bed`;
