#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Math::Complex;

# By Xilong CHEN
# Create date: 2022-08-22
# Contact: chen_xilong@outlook.com

# calculate the se of chr D

# write a summary of chr D, which is max and which is min

my $workdir =
  "/work/xchen/XCHEN_apple_genome/analysis_2022_Mar/run_Dsuite_genome";
my $input  = "$workdir/input";
my $output = "$workdir/output";
my $bin    = "$workdir/bin";
my $chr_f = "$input/s02.chr_list.txt";

open my $chr_r, "<", $chr_f or die "Could not open $chr_f:$!\n";
my @chrs;
while (my $line = <$chr_r>) {
    chomp $line;
    push @chrs, $line;
}
close $chr_f;

my $trio_f = "$input/Trios.txt";
my %trio_hash;
open my $in, "<", "$trio_f";
while (<$in>) {
    next if not /DOM/;
    chomp;
    $_ =~ s/\t/-/g;
    $trio_hash{$_} = 1;
}

# print Dumper \%trio_hash;

my $shdir = "s04.chr_D_summary";

# only use CID and DES trios
# plot for figure2 C

open my $ou1, ">", "$output/$shdir/chr_D_se_Rdata.txt";
open my $ou2, ">", "$output/$shdir/summary.txt";

print $ou1 "trio\tP2\tP3\tChr\tDstatistic\tSE\n";
my %chr_D_hash;

for my $chr (@chrs) {
    my $Dtree_f = "$output/s02.run_Dtest_by_chr/$chr" . "_tree.txt";
    open my $in_d, "<", $Dtree_f;
    <$in_d>;
    while (<$in_d>) {
        chomp;
        my @arr  = split /\t/, $_;
        my $trio = "$arr[0]-$arr[1]-$arr[2]";
        next if not exists $trio_hash{$trio};
        my $se = $arr[3] / $arr[4];
        print $ou1 "$trio\t$arr[1]\t$arr[2]\t$chr\t$arr[3]\t$se\n";
        $chr_D_hash{$trio}{$chr} = $arr[3];
    }
}

# print Dumper \%chr_D_hash;

open my $ou_g, ">", "$output/$shdir/Genome_D_se_Rdata.txt";
print $ou_g "trio\tP2\tP3\tChr\tDstatistic\tSE\n";
my %genome_hash;
open my $in_w, "<", "$output/s03.combine/genome_combined_tree.txt";
while (<$in_w>) {
    next if $. == 1;
    chomp;
    my @arr  = split /\t/, $_;
    my $trio = "$arr[0]-$arr[1]-$arr[2]";
    $genome_hash{$trio} = "$arr[3]\t$arr[4]";
    next if not exists $trio_hash{$trio};
    my $se = $arr[3] / $arr[4];
    print $ou_g "$trio\t$arr[1]\t$arr[2]\tGenome\t$arr[3]\t$se\n";
}

print $ou2 "trio\tgenome_D\tgenome_Z\tmax_chr\tmax_D\tmin_chr\tmin_D\n";
for my $trio ( sort keys %chr_D_hash ) {
    my @chrs =
      sort { $chr_D_hash{$trio}{$a} <=> $chr_D_hash{$trio}{$b} }
      keys %{ $chr_D_hash{$trio} };
    my $min_chr  = $chrs[0];
    my $max_chr  = $chrs[-1];
    my $min_D    = $chr_D_hash{$trio}{$min_chr};
    my $max_D    = $chr_D_hash{$trio}{$max_chr};
    my $genomeDZ = $genome_hash{$trio};

    print $ou2 "$trio\t$genomeDZ\t$max_chr\t$max_D\t$min_chr\t$min_D\n";
}
