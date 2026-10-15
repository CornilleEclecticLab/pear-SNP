#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create date: 2024-10-28
# Contact: chen_xilong@outlook.com

# ./analysis/fsc_modeling.demo_files/ApplePG1027/Group_G2/ApplePG_G2D2S2.run08/ApplePG_G2D2S2/ApplePG_G2D2S2.bestlhoods

# print "
# All my parameters have been set to output, so I can directly count the
# number of parameters and calculate AIC in bestlhoods files.
# If you hidden some parameters, you need to obtain the number of
# parameters from the *.est file.\n";

my $bestlhoods_path = shift @ARGV;
my $outpath         = shift @ARGV;
my @files           = glob($bestlhoods_path);

my %hash;
my @par_head;
for my $file (@files) {

    # print "$file\n";
    my $parent_dir           = dirname($file);
    my $grandparent_dir      = dirname($parent_dir);
    my $gg_dir               = dirname($grandparent_dir);
    my $bs_name              = basename($gg_dir);
    my $grandparent_dir_name = basename($grandparent_dir);
    my ( $model, $run ) = split /\./, $grandparent_dir_name;

    # $model =~ s/ApplePG_//g;
    open my $in, "<", $file;
    my $head = <$in>;
    chomp $head;
    @par_head = split /\t/, $head;
    pop @par_head;
    pop @par_head;
    # print  @par_head;
    my $line = <$in>;
    chomp $line;
    my @arr         = split /\t/, $line;
    my $MaxObsLhood = pop @arr;
    my $MaxEstLhood = pop @arr;
    my $n_par       = scalar @arr;

    # Calculate AIC
    my $log10_e = log( exp(1) ) / log(10);
    my $AIC     = 2 * $n_par - 2 * ( $MaxEstLhood / $log10_e );

    # Calculate delta-likelihood
    my $deltaL = $MaxObsLhood - $MaxEstLhood;

    $hash{$bs_name}{$run}{'MaxEstLhood'} = $MaxEstLhood;
    $hash{$bs_name}{$run}{'MaxObsLhood'} = $MaxObsLhood;
    $hash{$bs_name}{$run}{'AIC'}         = $AIC;
    $hash{$bs_name}{$run}{'deltaL'}      = $deltaL;
    $hash{$bs_name}{$run}{'par'}         = \@arr;
}
# print  @par_head;
# open my $ou1, ">", "$outpath/MaxEstLhood_AIC_boxplot.txt";

# print $ou1 "Boot_name\trun_nu\tMaxEstLhood\tAIC\n";
# for my $model ( sort keys %hash ) {
#     for my $run ( sort keys %{ $hash{$model} } ) {
#         print $ou1 "$model\t$run\t$hash{$model}{$run}{'MaxEstLhood'}\t$hash{$model}{$run}{'AIC'}\n";
#     }
# }

# find the best one for each model
my $par_head = join "\t", @par_head;

my %best_hash;
open my $ou2, ">", "$outpath/Best_run_for_each_boot.txt";
print $ou2 "Boot_name\trun_nu\tMaxEstLhood\tAIC\tdeltaL\t$par_head\n";
for my $model ( sort keys %hash ) {
    my @runs = sort {
        $hash{$model}{$b}{'MaxEstLhood'} <=> $hash{$model}{$a}{'MaxEstLhood'}
      }
      keys %{ $hash{$model} };
    print $ou2 "$model\t$runs[0]\t$hash{$model}{$runs[0]}{'MaxEstLhood'}\t";
    print $ou2
      "$hash{$model}{$runs[0]}{'AIC'}\t$hash{$model}{$runs[0]}{'deltaL'}\t";
    print $ou2 join "\t", @{ $hash{$model}{ $runs[0] }{'par'} };
    print $ou2 "\n";

    $best_hash{$model}{'MaxEstLhood'} =
      $hash{$model}{ $runs[0] }{'MaxEstLhood'};
    $best_hash{$model}{'MaxObsLhood'} =
      $hash{$model}{ $runs[0] }{'MaxObsLhood'};
    $best_hash{$model}{'AIC'}    = $hash{$model}{ $runs[0] }{'AIC'};
    $best_hash{$model}{'deltaL'} = $hash{$model}{ $runs[0] }{'deltaL'};
    $best_hash{$model}{'run'}    = $runs[0];
}

# # find the best of all based on MaxExtLhood and AIC
# open my $ou3, ">", "$outpath/Best_run_for_all.summary.txt";
# print $ou3 "Scenario\trun_nu\tMaxEstLhood\tAIC\tdeltaL\tbest_of_what\n";

# my @best_models_MaxEstLhood =
#   sort { $best_hash{$b}{'MaxEstLhood'} <=> $best_hash{$a}{'MaxEstLhood'} }
#   keys %best_hash;

# print $ou3 "$best_models_MaxEstLhood[0]\t";
# print $ou3 "$best_hash{$best_models_MaxEstLhood[0]}{'run'}\t";
# print $ou3 "$best_hash{$best_models_MaxEstLhood[0]}{'MaxEstLhood'}\t";
# print $ou3 "$best_hash{$best_models_MaxEstLhood[0]}{'AIC'}\t";
# print $ou3 "$best_hash{$best_models_MaxEstLhood[0]}{'deltaL'}\t";
# print $ou3 "best_of_MaxEstLhood\n";

# my @best_models_AIC =
#   sort { $best_hash{$a}{'AIC'} <=> $best_hash{$b}{'AIC'} }
#   keys %best_hash;

# print $ou3 "$best_models_AIC[0]\t";
# print $ou3 "$best_hash{$best_models_AIC[0]}{'run'}\t";
# print $ou3 "$best_hash{$best_models_AIC[0]}{'MaxEstLhood'}\t";
# print $ou3 "$best_hash{$best_models_AIC[0]}{'AIC'}\t";
# print $ou3 "$best_hash{$best_models_AIC[0]}{'deltaL'}\t";
# print $ou3 "best_of_AIC\n";

# # best par confidence interval data collection

# if ( $best_models_MaxEstLhood[0] eq $best_models_AIC[0] ) {
#     open my $ou4, ">", "$outpath/Best_MaxEstLhood_AIC_runs_pars.txt";
#     my $header_printed = 0;
#     for my $file (@files) {
#         my $parent_dir           = dirname($file);
#         my $grandparent_dir      = dirname($parent_dir);
#         my $grandparent_dir_name = basename($grandparent_dir);
#         my ( $model, $run ) = split /\./, $grandparent_dir_name;
#         $model =~ s/ApplePG_//g;
#         next if not $best_models_MaxEstLhood[0] eq $model;
#         open my $in_fh, "<", $file;
#         my $header = <$in_fh>;

#         if ( !$header_printed ) {
#             print $ou4 $header;
#             $header_printed = 1;
#         }
#         while ( my $line = <$in_fh> ) {
#             print $ou4 $line;
#         }
#         close $in_fh;
#         # if ( $best_hash{ $best_models_MaxEstLhood[0] }{'run'} eq $run ) {
#         #     `cp $file $outpath/Best_MaxEstLhood_AIC.bestlhoods`;
#         # }
#     }
# }
# else {
#     open my $ou5, ">", "$outpath/Best_MaxEstLhood_runs_pars.txt";
#     my $header_printed = 0;
#     for my $file (@files) {
#         my $parent_dir           = dirname($file);
#         my $grandparent_dir      = dirname($parent_dir);
#         my $grandparent_dir_name = basename($grandparent_dir);
#         my ( $model, $run ) = split /\./, $grandparent_dir_name;
#         $model =~ s/ApplePG_//g;
#         next if not $best_models_MaxEstLhood[0] eq $model;
#         open my $in_fh, "<", $file;
#         my $header = <$in_fh>;

#         if ( !$header_printed ) {
#             print $ou5 $header;
#             $header_printed = 1;
#         }
#         while ( my $line = <$in_fh> ) {
#             print $ou5 $line;
#         }
#         close $in_fh;
#         # if ( $best_hash{ $best_models_MaxEstLhood[0] }{'run'} eq $run ) {
#         #     `cp $file $outpath/Best_MaxEstLhood.bestlhoods`;
#         # }

#     }
#     open my $ou6, ">", "$outpath/Best_AIC_best_pars.txt";
#     $header_printed = 0;
#     for my $file (@files) {
#         my $parent_dir           = dirname($file);
#         my $grandparent_dir      = dirname($parent_dir);
#         my $grandparent_dir_name = basename($grandparent_dir);
#         my ( $model, $run ) = split /\./, $grandparent_dir_name;
#         $model =~ s/ApplePG_//g;
#         next if not $best_models_AIC[0] eq $model;
#         open my $in_fh, "<", $file;
#         my $header = <$in_fh>;

#         if ( !$header_printed ) {
#             print $ou6 $header;
#             $header_printed = 1;
#         }
#         while ( my $line = <$in_fh> ) {
#             print $ou6 $line;
#         }
#         close $in_fh;
#         # if ( $best_hash{ $best_models_AIC[0] }{'run'} eq $run ) {
#         #     `cp $file $outpath/Best_AIC.bestlhoods`;
#         # }
#     }
# }
