#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Temp qw/tempfile/;

# By Xilong CHEN
# Create date: 2025-06-10
# Contact: chen_xilong@outlook.com

my $file_bed    = shift;    # masked bed
my $file_tsv    = shift;    # tsv file
my $filt_treegz = shift;    # tree file with csv
my $cut_off     = shift;

# csv to bed file
my ( $ou_tb, $temp_bed ) =
  tempfile( SUFFIX => '.bed', DIR => '.', UNLINK => 1 );
open my $in_csv, "<", $file_tsv;
my %line_nu_hash;
while (<$in_csv>) {
    chomp;
    next if $. == 1;
    my ( $chr, $start, $end ) = split /\t/, $_;
    my $bs = $start - 1;
    my $be = $end;
    print $ou_tb "$chr\t$bs\t$be\n";
    my $key = join( "\t", $chr, $bs, $be );
    $line_nu_hash{$key} = $. - 1;
}
close $ou_tb;

open my $in, "-|", "bedtools intersect -a $file_bed -b $temp_bed -wo"
  or die "Failed to run bedtools: $!\n";
my %overlap;
my %length;
while (<$in>) {
    chomp;
    my @fields = split(/\t/);
    my ( $chr, $start, $end, $overlap_len ) = @fields[ -4, -3, -2, -1 ];
    my $key = join( "\t", $chr, $start, $end );
    $overlap{$key} += $overlap_len;
    $length{$key} = $end - $start;
}
close $in;
my %rm_hash;
foreach my $key ( sort keys %overlap ) {
    my $frac = $overlap{$key} / $length{$key};
    if ( $frac >= $cut_off ) {
        my $index = $line_nu_hash{$key};
        $rm_hash{$index} = 1;
    }
}

open my $in_tree, "-|", "gunzip -c $filt_treegz";
while (<$in_tree>) {
    if ( not exists $rm_hash{$.} ) {
        print $_;
    }
}

