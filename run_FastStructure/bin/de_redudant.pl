#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw(sum);

# By Xilong CHEN
# Create date: 2022-03-10
# Contact: chen_xilong@outlook.com

my @Qfiles = glob "./input_for_pop/*.Q";

my $dir = "de_redundant_qfile";
print "Creat folder $dir\n";
`rm -rf $dir`;
`mkdir $dir`;

# print join "\n", @Qfiles;
for my $f (@Qfiles) {
    my %q_hash;
    my %sum_hash;
    my $new_k_nu;
    open my $in1, "<", "$f";
    while (<$in1>) {
        chomp;
        my @arr = split / /, $_;
        for my $i ( 0 .. $#arr ) {
            my $k = $i + 1;
            push @{ $q_hash{$k} }, $arr[$i];
        }
    }
    for my $k ( sort { $a <=> $b } keys %q_hash ) {
        my $sum = sum @{ $q_hash{$k} };
        if ( $sum != 0 ) {
            $sum_hash{$k} = 1;
            $new_k_nu++
        }
    }
    open my $in2, "<", "$f";
    # print Dumper \%sum_hash;
    # print $new_k_nu;
    my $fname = basename $f;
    open my $ou, ">", "$dir/$fname.$new_k_nu.meanQ";
    while (<$in2>) {
        chomp;
        my @arr = split / /, $_;
        my @new_arr;
        for my $i ( 0 .. $#arr ) {
            my $k = $i + 1;
            push @new_arr, $arr[$i] if exists $sum_hash{$k};
        }
        print $ou join " ", @new_arr;
        print $ou "\n";
    }
    # print Dumper \%hash;
}
