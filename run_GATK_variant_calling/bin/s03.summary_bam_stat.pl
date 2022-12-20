#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/14/2020 16:06
# Contact: chen_xilong@outlook.com

my $workdir = "~/work/pear/run_GATK_variant_calling";
my @bamreads = glob "$workdir/output/s02.stat_bam_by_bbmap/*.bamreads";


my $shdir = "s03.summary_bam_stat";
my $output  = "$workdir/output/$shdir/bam.stat.txt";
open my $ou, ">", $output;

print $ou "ID\tMapped Reads (M)\tMapped bases\tmapped Depth (x)\tMapped Coverage (%)\n";
for (@bamreads) {
    my ($name) = fileparse ($_, qr/\.[^.]*/);
    my ($mr, $mb, $ac, $pc);
    open my $fh, "<", "$_";
    while (<$fh>) {
        chomp;
        if (/Mapped reads/) {
            my @line = split /\s+/, $_;
            $mr = $line[-1];
        } elsif (/Mapped bases/) {
                        my @line = split /\s+/, $_;
            $mb = $line[-1];
        } elsif (/Average coverage/) {
                        my @line = split /\s+/, $_;
            $ac = $line[-1];
        } elsif (/Percent of reference bases covered/) {
                        my @line = split /\s+/, $_;
            $pc = $line[-1];
        }
    }
    print $ou join "\t", ($name, $mr, $mb, $ac, $pc);
    print $ou "\n";
}




# Mapped reads:                           200291577
# Mapped bases:                           17415291175
# Average coverage:                       24.524
# Percent of reference bases covered:     86.42


# print join "\n", @stats;
