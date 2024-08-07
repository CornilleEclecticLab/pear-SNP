#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;

# By Chen Xilong
# Create Date: 2/6/2020 15:29
# Contact: chen_xilong@outlook.com

# Updated by Yuqi to get the double digit K in the output file name.
# Updated the regex by Yuqi to capture the correct name when 
#   clumpak result folder name ends with a digital.

my $usage = "perl $0 <the clumpak_result folder> <min K number> <max K number> \n\n";
die $usage if @ARGV != 3;

my $re_folder = shift @ARGV;
my $min = shift @ARGV;
my $max = shift @ARGV;

`mkdir input_for_pophelper`;
for my $k_f ( $min .. $max ) {
    my @Qfiles
        = glob
        "$re_folder/K=$k_f/M*Cluster*/CLUMPP.files/ClumppIndFile.output";
    # print join "\n", @Qfiles;
    for my $f (@Qfiles) {
        $f =~ /K=\d+\/(.+?Cluster\d*)/;
        my $m = $1;
        open my $in, "<", "$f";
        my $format_kf = sprintf "%02d", $k_f;
        print "$k_f/K.$format_kf.$m.Q\n";
        open my $ou, ">", "./input_for_pophelper/K.$format_kf.$m.Q";
        while (<$in>) {
            $_ =~ s/^.+: //;
            print $ou $_;
        }
    }

}

