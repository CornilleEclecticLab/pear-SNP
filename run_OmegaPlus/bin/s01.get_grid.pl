#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2021-10-18
# Contact: chen_xilong@outlook.com

my $usage = "perl $0 <window_size> <karyotype file>\n";
die $usage if $#ARGV != 1;

my @LOG = @ARGV;
my $window =  shift @ARGV;
my $karyo_file = shift @ARGV;

my $time = `date +"%y%m%d-%H%M"`;
chomp $time;
open my $oulog, ">", "$0.$time.log";
print $oulog "$usage\nperl $0 ";
print $oulog join ' ', @LOG;
print $oulog "\n";

my $workdir =
  "/shared/home/ynie/work/pear/run_RAiSD/";
my $input = "$workdir/input";

my %karyo_hash;
open my $kin, "<", $karyo_file or die;
open my $ou,         ">", "$input/grid_chr.$window.txt";
while (<$kin>) {
#    next if /Chr00/;
    my ( $chr, $l ) = split /\t/, $_;
    my $grid = int( $l / $window ) + 1;
    print $ou "$chr\t$grid\n";
}