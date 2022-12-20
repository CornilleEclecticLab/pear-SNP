#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/14/2020 11:02
# Contact: chen_xilong@outlook.com
# Adapted by Yuqi NIE, on 2022-10-27. Changed the memory used for GATK

my $shdir = "s04.combine_chr_gvcf_to_vcf_allsite";
`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;

my $bwa2     = "module load bioinfo/bwa-mem2-2.0";
my $samtools = "module load bioinfo/samtools-1.9";
my $gatk     = "module load bioinfo/gatk-4.1.7.0";
my $workdir  = "~/work/pear/run_GATK_variant_calling";

my $output_path = "$workdir/output/$shdir";
my $genome      = "$workdir/input/Genome_reference/Prunus_persica_v2.0.a1_scaffolds.fasta";

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
    open my $ou2, ">", "$chr.gvcf.list";
    my @files = glob "$workdir/output/s01.fastq_to_gvcf.peach/*/*$chr*g.vcf.gz";
    for (@files) {
        next if /A326|A352/;
    print $ou2 "$_\n";
    }

}

for my $chr (@Chrs) {
    my $file         = "$workdir/bin/$chr.gvcf.list";
    my $job          = "$chr.s04";
    my $combine_gvcf = <<"END_COM";
#!/bin/bash
#SBATCH -J $job
#SBATCH -o $job.out
#SBATCH -e $job.err
module purge
$samtools
$gatk

gatk --java-options "-Xmx32g" CombineGVCFs \\
    -R $genome \\
    -V $file \\
    -O $output_path/$chr.combine.g.vcf.gz && \\
echo "combine gvcfs done!"

gatk --java-options "-Xmx32g" GenotypeGVCFs \\
    -all-sites \\
    -R $genome \\
    -V $output_path/$chr.combine.g.vcf.gz \\
    -O $output_path/$chr.combine.vcf.gz && \\
echo "Genotype gvcf to vcf done!"

END_COM

    open my $ou, ">", "$shdir/$job.sh" or die;

    print $ou $combine_gvcf;

}

