#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use File::Basename;

# By Chen Xilong
# # Create Date: 2/11/2020 11:17
# # Contact: chen_xilong@outlook.com
# Adapted by Yuqi on 2022-10-24

open my $ou, ">", "peach.list";
print $ou join "\n", glob "/work/zruilin/Peach/run_QC/input/peach_Illumina_novogene_25082022/merged_row_data/*/*.fq.gz";
print $ou "\n";
print $ou join "\n", glob "/work/zruilin/Peach/run_QC/input/peach_Illumina_novogene_25082022/merged_row_data/*/*.fastq.gz";
print $ou "\n";
