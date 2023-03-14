#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw(shuffle);

# By Xilong CHEN
# Create date: 2021-12-20
# Contact: chen_xilong@outlook.com

# update 2023-03-14

my $usage = "perl $0 <window_size> <input_vcf.gz> <output_vcf>\n";

my $win     = shift @ARGV;
my $invcfgz = shift @ARGV;
my $ouf     = shift @ARGV;

open my $ou, ">", "$ouf";
open my $z, '-|', '/usr/bin/gunzip', '-c', "$invcfgz";

my @buffer;
while (<$z>) {
    if ( $_ =~ /^#/ ) {
        print $ou $_;
    }
    else {
        push @buffer, $_;
        if ( @buffer == $win || eof ) {
            @buffer = shuffle(@buffer);
            print $ou $buffer[0];
            @buffer = ();
        }
    }
}

