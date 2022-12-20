#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/14/2020 11:36
# Contact: chen_xilong@outlook.com

# not filter the DP
# updated 2022 Mar
# first filter all the site don't have QUAL, FILTER, INFO

my $shdir = "s05.hard_filter_chr_vcf_variant_invarinat";
`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;

my $bwa2     = "module load bioinfo/bwa-mem2-2.0";
my $samtools = "module load bioinfo/samtools-1.9";
my $gatk     = "module load bioinfo/gatk-4.1.7.0";
my $bcftools = "module load bioinfo/bcftools-1.9";
my $workdir  = "/work/zruilin/Peach/run_GATK_variant_calling";

my $output_04 = "$workdir/output/s04.combine_chr_gvcf_to_vcf_allsite";
my $output_path = "$workdir/output/$shdir";
my $genome_path = "$workdir/input/Genome_reference/Prunus_persica_v2.0.a1_scaffolds.fasta";

my @Chrs = qw(
Pp01
Pp02
Pp03
Pp04
Pp05
Pp06
Pp07
Pp08
scaffolds.list);


for my $chr (@Chrs) {
    my $job = "$chr.s05";
    my $filter = <<"END_FILTER";
#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err
module purge
$samtools
$gatk
$bcftools
# filter all the site don't have QUAL, FILTER, INFO by filter F_MISSING = 1

bcftools filter -e 'F_MISSING = 1' \\
$output_04/$chr.combine.vcf.gz -O z \\
-o $output_path/$chr.rmQFI.combine.vcf.gz

tabix $output_path/$chr.rmQFI.combine.vcf.gz

# you can use -XL to remove scaffold like mitochondria
gatk --java-options "-Xmx4g" VariantFiltration \\
    -R $genome_path \\
    -V $output_path/$chr.rmQFI.combine.vcf.gz \\
    -filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0) " \\
    --filter-name "badSeq" \\
    --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \\
    --filter-name "badStrand" \\
    --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \\
    --filter-name "badMap" \\
    -O $output_path/$chr.filted_anno.vcf.gz && \\
echo "hard filter done!"

gatk --java-options "-Xmx4g" SelectVariants \\
-R $genome_path \\
-V $output_path/$chr.filted_anno.vcf.gz \\
--exclude-filtered true \\
-O $output_path/$chr.filted_passed_sites.vcf.gz && \\
echo "generate the passed vcf done!"


END_FILTER
    open my $ou1, ">", "$shdir/$job.sh" or die;
    print $ou1 $filter;
}
