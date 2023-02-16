#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 12/16/2020 22:44
# Contact: chen_xilong@outlook.com

my $usage = "perl $0 <nex file> <color information file> <output file>
The script will color the lable and the closted edge for splitstree4 nex file.
1 open your tree by splitstree.
2 save the file as nex file.
3 prepare color information file.
4 run the script.
5 open the output nex file by splitstree.
6 print it as a PDF file.
The <color information file> have three conlums(ID  group   color). Separator is TAB. eg,
CULT_ROOT_DOM_USA_ncA00<tab>CULT_DOM<tab>#E41A1C

";

die $usage if @ARGV != 3;
my ( $nex_file, $color_file, $out_file ) = @ARGV;
# my $nex_file   = "test.dist2.nex";
# my $color_file = "group_color.txt";
# my $out_file   = "test.dist2.color.nex";

open my $in1, "<", "$color_file";
my %c_hash;
while (<$in1>) {
    chomp;
    my ( $id, $gourp, $color ) = split /\t/, $_;
    # $color =~ s/#//;
    $c_hash{$id} = $color;
}



my %line_nu;

# find the VLABELS line number

open my $in_n1, "<", "$nex_file";
while (<$in_n1>) {
    if ( $_ =~ /^VLABELS/ ) {
        $line_nu{'VLABELS'} = $.;
    }
    if ( exists $line_nu{'VLABELS'} and $_ =~ /;/ ) {
        $line_nu{'VLABELS_END'} = $.;
        last;
    }
}

open my $in_n2, "<", "$nex_file";
while (<$in_n2>) {
    if ( $_ =~ /^EDGES/ ) {
        $line_nu{'EDGES'} = $.;
    }
    if ( exists $line_nu{'EDGES'} and $_ =~ /;/ ) {
        $line_nu{'EDGES_END'} = $.;
        last;
    }
}

print Dumper \%line_nu;

# build the id_nu => id hash
my %color_hash;

open my $in_r1, "<", "$nex_file";
while (<$in_r1>) {
    if ( $. > $line_nu{'VLABELS'} and $. < $line_nu{'VLABELS_END'} ) {
        chomp;
        my @arr = split / /, $_;
        my $id  = $arr[0];
        $arr[1] =~ s/'//g;
        my $name = $arr[1];

        # print "$1\t$2\n";
        $color_hash{$id} = $c_hash{$name};
    }

}

# print Dumper \%color_hash;
# exit;
open my $ou,    ">", "$out_file";
open my $in_r2, "<", "$nex_file";
while (<$in_r2>) {
    if ( $. > $line_nu{'VLABELS'} and $. < $line_nu{'VLABELS_END'} ) {
        chomp;
        my @arr = split / /, $_;
        my $id  = $arr[0];
        $_ =~ s/,/ lc=$color_hash{$id},/;
        print $ou "$_\n";
    }
    elsif ( $. > $line_nu{'EDGES'} and $. < $line_nu{'EDGES_END'} ) {
        chomp;
        my @arr = split / /, $_;
        my $id  = $arr[2];
        if ( exists $color_hash{$id} ) {
            $_ =~ s/,/ fg=$color_hash{$id},/;
            print $ou "$_\n";
        }
        else {
            print $ou "$_\n";
        }
    }
    else {
        print $ou $_;
    }
}

# 235 'cropIR_wild_hyb' x=-85 y=4,
# 235 'cropIR_wild_hyb' x=-85 y=4 lc=#00FFFF,

# close $in;

# open my $in2, "<", "$nex_file";
# while (<$in2>) {

# }

# my @keys   = keys %c_hash;

# my $nu_cut = @keys + 1;

# open my $in, "<", "$nex_file";
# open my $ou, ">", "$out_file";

# my %nu_hash;
# while (<$in>) {
#     if ( $_ =~ /^\tN:/ ) {
#         $_ =~ /^\tN: (\d+) /;
#         my $nu = $1;
#         if ( $nu < 2 or $nu > $nu_cut ) {
#             print $ou $_;
#         }
#         else {
#             $_ =~ /(0x[\w]{6}).+: '(.+)' .+(0x[\w]{6}).+/;
#             my $node_color  = $1;
#             my $lable_color = $3;
#             my $id          = $2;
#             my $color       = $c_hash{$id};
#             $nu_hash{$nu} = $color;
#             $_ =~ s/0x[\w]{6}/0x$color/g;
#             print $ou $_;
#         }

#     }
#     elsif ( $_ =~ /^\tE:/ ) {
#         $_ =~ /^\tE: (\d+) (\d+) /;
#         my $nu = $2;
#         if ( $nu < 2 or $nu > $nu_cut ) {
#             print $ou $_;
#         }
#         else {
#             $_ =~ /(0x[\w]{6})/;
#             my $edge_color = $1;
#             my $color      = $nu_hash{$nu};
#             $_ =~ s/0x[\w]{6}/0x$color/g;
#             print $ou $_;
#         }

#     }
#     else {
#         print $ou $_;
#     }

# }
