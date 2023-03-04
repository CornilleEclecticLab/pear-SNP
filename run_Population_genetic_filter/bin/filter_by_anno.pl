#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 10/22/2020 17:05
# Contact: chen_xilong@outlook.com

my $usage = "perl $0 <input.vcf.gz> <output.vcf>\n
Please bgzip and tabix the output vcf file\n";

die $usage if @ARGV != 2;

my $invcfgz = shift @ARGV;
my $ouf = shift @ARGV;

open my $ou, ">", "$ouf";
open my $z, '-|', '/usr/bin/gunzip', '-c', "$invcfgz";

while ( my $line = <$z> ) {
    if ( $line =~ /^#/ ) {
            print $ou $line;
                }
                    else {
                            if ( $line =~ /synonymous_variant/ ) {
                                        print $ou $line;
                                                }
                                                    }
                                                    }

                                                    close $ou;


