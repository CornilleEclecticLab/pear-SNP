#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;

# By Chen Xilong
# Create Date: 2/6/2020 15:29
# Contact: chen_xilong@outlook.com
my $usage = "perl $0 <the clumpak_result folder> <min K number> <max K number> \n\n";
die $usage if @ARGV != 3;

my $re_folder = shift @ARGV;
my $min =shift @ARGV;
my $max = shift @ARGV;

`mkdir input_for_pop`;
for my $k_f ( $min .. $max ) {
    my @Qfiles
        = glob
        "$re_folder/K=$k_f/M*Cluster*/CLUMPP.files/ClumppIndFile.output";
    # print join "\n", @Qfiles;
    for my $f (@Qfiles) {
        $f =~ /\d\/(.+?Cluster\d*)/;
        my $m = $1;
        open my $in, "<", "$f";
        print "$k_f/K.$k_f.$m.Q\n";
        open my $ou, ">", "./input_for_pop/K.$k_f.$m.Q";
        while (<$in>) {
            $_ =~ s/^.+: //;
            print $ou $_;
        }
    }

}

