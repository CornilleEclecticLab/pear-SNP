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

$/ = "\/\/\n";
open my $in, "<", "$f";
my $line = <$in>;
chomp $line;
print $line;
while (<$in>) {
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
    my @newp = sort { $a <=> $b } keys %hash;

    my $new_sigs = $#newp + 1;
    print "\n\/\/\n";
    print "segsites: $new_sigs\n";
    print "positions: ";
    print join " ", @newp;
    print "\n";
    for my $line (@block) {
        print substr($line, 0, $new_sigs);
        print "\n";
    } 
}

