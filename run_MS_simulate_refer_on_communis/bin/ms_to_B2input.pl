#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-05-25
# Contact: chen_xilong@outlook.com

my $f = shift @ARGV;
my $outpath = shift @ARGV;

my $name = basename ($f, ".txt");

$/ = "\/\/\n";
open my $in, "<", "$f";
my $line = <$in>;
chomp $line;

my $cm_per_bp = 1.92876E-06;


# print $line;
while (<$in>) {
    # next if $. > 3;
    chomp;
    my @block = split /\n/, $_;
    my $sigs  = shift @block;
    my $posi  = shift @block;
    my @arr   = split /\s+/, $posi;
    shift @arr;
    my %hash;

    for my $p (@arr) {
        $hash{$p} = 1;
    }

    # print $posi;
    my @newp      = sort { $a <=> $b } keys %hash;
    my $end_index = length $block[0] - 1;

    my @matrix;
    for my $line (@block) {
        my @arr = split //, $line;
        push @matrix, \@arr;
    }

    my @tr;
    for my $row ( 0 .. @matrix - 1 ) {
        for my $col ( 0 .. @{ $matrix[$row] } - 1 ) {
            $tr[$col][$row] = $matrix[$row][$col];
        }
    }
    my $re = $. - 1;
    open my $ou, ">", "$outpath/B2input.$name.re$re.txt";
    for my $i ( 0 .. $#newp ) {
        my $po  = 1000000 * $newp[$i];
        my $str = join "", @{ $tr[$i] };
        my $n   = length $str;
        my $x   = ( $str =~ tr/1/1/ );
        my $gd = $po * $cm_per_bp;
        print $ou "$po\t$gd\t$x\t$n\n";
    }
}
