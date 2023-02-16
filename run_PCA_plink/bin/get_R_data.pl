#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 3/29/2021 10:31
# Contact: chen_xilong@outlook.com

open my $in, "<", "Group_subset.txt" or die;
my %g_hash;

while (<$in>) {
    chomp;
    my ($id, $g) = split /\t/, $_;
    $g_hash{$id} = $g;
}
close $in;

open my $ou, ">", "pca_input.data";
my @head;
push @head, "Group";
for my $i (1..20) {
    push @head, "PC$i";
}
print $ou join "\t", @head;
print $ou "\n";

open my $in2, "<", "plink.eigenvec" or die;
while (<$in2>) {
    chomp;
    my @line = split / /, $_;
    my $g = $g_hash{$line[1]};
    print $ou "$g\t";
    print $ou join "\t", @line[2..21];
    print $ou "\n";
}

