#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 4/6/2021 10:23
# Contact: chen_xilong@outlook.com

`mkdir -p ../output/s03.halfR2_and02R2_keepRbin`;

my $input = "../output/s01.generate_PopLDdecay_running_scripts";
open my $ou, ">", "../output/s03.halfR2_and02R2_keepRbin/result.txt";

my @file_list = glob ("$input/*.bin");

print $ou Dumper \@file_list;

print $ou "\n";
print $ou "pop\tlength_when_r2=1/2maxr2\n";
for my $f (@file_list) {
    # next if $f =~ /All/;
    open my $fh_bgz_in, "<", "$f";
    <$fh_bgz_in>;
    my $line1 = <$fh_bgz_in>;
    chomp $line1;
    my ( $l1, $r_sq1 ) = split /\t/, $line1, 3;
    while ( my $line = <$fh_bgz_in> ) {
        chomp $line;
        my ( $l, $r_sq ) = split /\t/, $line, 3;
        if ( ( $r_sq - ( $r_sq1 / 2 ) ) < 0 ) {
            print $ou "$f\t$l\n";
            last;
        }
    }
}

print $ou "\n";
print $ou "pop\tlength_when_r2<0.2\n";
for my $f (@file_list) {
    # next if $f =~ /All/;
    open my $fh_bgz_in, "<", "$f";
    <$fh_bgz_in>;
    my $line1 = <$fh_bgz_in>;
    chomp $line1;
    my ( $l1, $r_sq1 ) = split /\t/, $line1, 3;
    while ( my $line = <$fh_bgz_in> ) {
        chomp $line;
        my ( $l, $r_sq ) = split /\t/, $line, 3;
        if ( $r_sq < 0.2 ) {
            print $ou "$f\t$l\n";
            last;
        }
    }
}

