#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-11-29
# Contact: chen_xilong@outlook.com

my $set = shift @ARGV;

my $rel_file      = "$set.pca.rel";
my $eigenval_file = "$set.pca.eigenval";
my $pve_file      = "$set.pca.pve";

my $sum;
open my $in_rel, "<", $rel_file;
while (<$in_rel>) {
    chomp;
        my @arr = split /\t/, $_;
            $sum += $arr[-1];
            }

            open my $in_eigenval, "<", $eigenval_file;
            open my $ou_pve,      ">", $pve_file;
            while (<$in_eigenval>) {
                chomp;
                    my $l_pve = $_ / $sum * 100;
                        print $ou_pve "$l_pve\n";
                        }


