#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-03-31
# Contact: chen_xilong@outlook.com

my $workdir = "/work/zruilin/Peach/QC";
my $output  = "$workdir/output";
my $input   = "$workdir/input";

my $usage = "perl $0 <set_list_file>
change fastp optins by differnt set raw qc report\n";
die $usage if @ARGV != 1;

my $set = shift @ARGV;
my $set_name = basename $set;

my $shdir = "s03.fastp_and_qc.$set_name";
print "creat $shdir
use 2 cpu\n";

`rm -rf $shdir`;
`mkdir $shdir`;

open my $in, "<", "$set" or die;
my %Dnu;
while (<$in>) {
    chomp;
    my @aa  = fileparse $_;
    my $dnu = basename $aa[1];
    # print "$dnu\n";
    push @{ $Dnu{$dnu} }, $_;
}

# print Dumper \%Dnu;

for my $d ( sort keys %Dnu ) {
    # print "$d\n";
    my $job = "$d.$shdir";
    open my $ou, ">", "$shdir/$job.sh";
    print $ou "#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err
module purge
  
module load bioinfo/FastQC_v0.11.7
module load bioinfo/fastp-0.21.0

mkdir $output/$shdir/$d
";
    my @f = sort @{ $Dnu{$d} };
    my ( @R1, @R2 );
    for my $i ( 0 .. $#f ) {
        if ( $i % 2 eq 0 ) {
            push @R1, $f[$i];
        }
        else {
            push @R2, $f[$i];
        }
    }
    for my $j ( 0 .. $#R1 ) {
        my $in_f1 = $R1[$j];
        my $in_f2 = $R2[$j];
        print "$in_f1\t$in_f2\n";
        my $f1    = basename($in_f1);
        my $f2    = basename($in_f2);
        my $sh    = << "SH_END";
fastp -f 2 -l 50 \\
-i $in_f1 -I $in_f2 \\
-o $output/$shdir/$d/clean.$f1 -O $output/$shdir/$d/clean.$f2 \\
-h $output/$shdir/$d/fastp.$d.html -j $output/$shdir/$d/fastp.$d.json -R '$d' &&\\
fastqc -t 2 -f fastq -o $output/$shdir/$d $output/$shdir/$d/clean.$f1 &&\\
fastqc -t 2 -f fastq -o $output/$shdir/$d $output/$shdir/$d/clean.$f2 &&\\
echo done

SH_END
        print $ou $sh;
    }

}
