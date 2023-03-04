#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;

# By Chen Xilong
# Create Date: 3/27/2018 11:28
# Contact: chen_xilong@outlook.com

my @files = glob './log/*.log';

# print Dumper \@files;

my %cv;
for my $f (@files) {
    $f =~ /\.(\d+)\./;
    my $K = $1;
    open my $in, "<", "$f" or die;
    while (<$in>) {
        chomp;
        next if not /CV error/;
        $_ =~ s/CV error = //;
        my ( $cv_error, $un ) = split /\,/, $_;
        push @{ $cv{$K} }, $cv_error;
    }
}

# print Dumper \%cv;
open my $out, ">", 'fastst_cv.data' or die;
print $out "K\tCV_error\n";
for my $K ( sort { $a <=> $b } keys %cv ) {
    for my $c ( @{ $cv{$K} } ) {
        print $out "$K\t$c\n";
    }

}
