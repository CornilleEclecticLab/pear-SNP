#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2022-08-02
# Contact: chen_xilong@outlook.com

# Adjusted by Yuqi
# 2024-09-04 

my $folder = "s03.find_FPRcutoff_OmegaPlus";
`mkdir -p ../output/$folder`;

my $bindir = dirname(abs_path($0));
my $workdir = dirname($bindir);
my $input  = "$workdir/input";
my $output = "$workdir/output";

my $s02 = "$workdir/output/s02.generate_positive_selection_detect_OmegaPlus";


open my $in, "<", "$input/s01.sample_tab_population.txt" or die;
my %pop;
while (<$in>) {
    chomp;
    my ( $id, $group ) = split /\t/, $_;
    $pop{$group} += 1;
}
my @pops = keys %pop;

my $cutoff    = 0.05;      # FPR at 95%
my $cutoff_op = 1 - $cutoff;




for my $pop (@pops) {
    my @omegas;
    for my $re ( 1 .. 10 ) {
        # $/ = "\n\/\/";

        # open my $inSw, "<", "$s02/SweeD_Report.$pop.$re";
        # <$inSw>;
        # while (<$inSw>) {
        #     chomp;
        #     my @block = split /\n/, $_;
        #     my $nu    = shift @block;
        #     my $head    = shift @block;
        #     # print "$nu\n";
        #     my @clr_re;
        #     for my $line (@block) {
        #         my ( $po, $clr ) = split /\t/, $line;
        #         push @clr_re, $clr;
        #     }
        #     @clr_re = reverse sort { $a <=> $b } @clr_re;
        #     push @clrs, $clr_re[0];
        # }
        $/ = "\n\/\/";
        open my $inOm, "<", "$s02/OmegaPlus_Report.$pop.REP$re";
        <$inOm>;
        while (<$inOm>) {
            chomp;
            my @block = split /\n/, $_;
            my $nu    = shift @block;
            # my $head    = shift @block;
            # print "$nu\n";
            my @omega_re;
            for my $line (@block) {
                my ( $po, $omega ) = split /\t/, $line;
                push @omega_re, $omega;
            }
            @omega_re = reverse sort { $a <=> $b } @omega_re;
            push @omegas, $omega_re[0];
        }

        # $/ = "\/\/";
        # open my $inRA, "<", "$s02/RAiSD_Report.$pop.$re";
        # <$inRA>;
        # while (<$inRA>) {
        #     chomp;
        #     my @block = split /\n/, $_;
        #     my $nu    = shift @block;
        #     # print "$nu\n";
        #     my @mu_re;
        #     for my $line (@block) {
        #         my ( $po, $mu ) = split /\t/, $line;
        #         push @mu_re, $mu;
        #     }
        #     @mu_re = reverse sort { $a <=> $b } @mu_re;
        #     push @mus, $mu_re[0];
        # }

    }

    # @clrs = reverse sort { $a <=> $b } @clrs;
    # pop @clrs for ( 1 .. int( @clrs * $cutoff_op ) );
    # my $clrs_cutoff = $clrs[-1];

    # # biger -> positive
    @omegas = reverse sort { $a <=> $b } @omegas;
    my $count = scalar @omegas;
    print "The number of elements in \@omegas is $count.\n";
    pop @omegas for ( 1 .. int( @omegas * $cutoff_op ) - 1 );
    my $omegas_cutoff = $omegas[-1];
    print "$pop\t$omegas_cutoff\n";

    # biger -> positive
    # @mus = reverse sort { $a <=> $b } @mus;
    # # print Dumper \@mus;
    # pop @mus for ( 1 .. int( @mus * $cutoff_op ) );
    # my $mus_cutoff = $mus[-1];
    # print "$pop\t$mus_cutoff\n";
    open my $ouc, ">", "$output/$folder/$pop.FPRcutoff.Rdata.txt";
    print $ouc "omegas_cutoff $pop\n";
    print $ouc "$omegas_cutoff\n";
}

