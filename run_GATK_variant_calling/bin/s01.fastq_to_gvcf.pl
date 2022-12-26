#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

# By Xilong CHEN
# Create Date: 8/7/2020 01:02
# Contact: chen_xilong@outlook.com

my $step = "s01.fastq_to_gvcf";
my $group = "pear";
my $shdir = "$step.$group";

`rm -rf $shdir` if ( -d "$shdir" );
`mkdir $shdir`;
my $workdir  = "~/work/pear/run_GATK_variant_calling/";

open my $sou, ">", "$shdir/00sbatch.sh";

my $filelist = "$group.list";

open my $in, "<", $filelist or die;

my %Dnu;
while (<$in>) {
    chomp;
    my @aa  = fileparse $_;
    my $dnu = basename $aa[1];
    print "$dnu\n";
    push @{ $Dnu{$dnu} }, $_;
}

my $bwa2     = "module load bioinfo/bwa-mem2-2.0";
my $samtools = "module load bioinfo/samtools-1.9";
my $gatk     = "module load bioinfo/gatk-4.1.7.0";

my @Chrs = qw(
GWHBAOS00000076
GWHBAOS00000158
GWHBAOS00000386
GWHBAOS00000381
GWHBAOS00000224
GWHBAOS00000085
GWHBAOS00000163
GWHBAOS00000172
GWHBAOS00000128
GWHBAOS00000425
GWHBAOS00000352
GWHBAOS00000424
GWHBAOS00000274
GWHBAOS00000099
GWHBAOS00000356
GWHBAOS00000365
GWHBAOS00000335
"$workdir/input/input/scaffolds.list");


my $out    = "$workdir/output/$shdir";
my $genome = "$workdir/input/Genome_reference/Prunus_persica_v2.0.a1_scaffolds.fasta";

for my $d ( sort keys %Dnu ) {
    print "$d\n";
    open my $ou, ">", "$shdir/$d.mapping.sh" or die;
    my $jobmap = $d."map";
    print $ou "#!/bin/bash\n";
    print $ou
        "if [ -d \"$out/$d\" ]; then rm -rf  $out/$d; fi;mkdir $out/$d\n";
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
    my $R1_list = join " ", @R1;
    my $R2_list = join " ", @R2;
    my $mapping = <<"END_MAPPING";

#SBATCH -J $jobmap
#SBATCH -o $jobmap.out
#SBATCH -e $jobmap.err

module purge

$bwa2
$samtools

bwa-mem2.sse41 mem -t 4 $genome \\
'<zcat $R1_list' \\
'<zcat $R2_list' \\
-R '\@RG\\tID:$d\\tSM:$d' \\
| samtools view -S -b - > $out/$d/$d.bam

END_MAPPING

    print $ou "$mapping\n";
    print $sou "submap$d=\$(sbatch -c 4 --mem=20G $d.mapping.sh)\njidm$d=\$(echo \$submap$d | awk -F' ' '{print \$4}')\n";

    open my $ou2, ">", "$shdir/$d.sort_markdu_bam.sh" or die;
    my $jobbam = $d."bam";
    my $bam = <<"END_BAM";
#!/bin/bash
#SBATCH -J $jobbam
#SBATCH -o $jobbam.out
#SBATCH -e $jobbam.err
module purge
$samtools
$gatk
samtools sort -@ 4 -m 2G \\
-O bam -o $out/$d/$d.sorted.bam $out/$d/$d.bam && echo "sort bam done!"

rm $out/$d/$d.bam

gatk --java-options "-Xmx8g" MarkDuplicates -I $out/$d/$d.sorted.bam \\
-O $out/$d/$d.sorted.markdu.bam \\
-M $out/$d/$d.sorted.markdu_metrics.txt \\
&& echo "Mark dupilcate done!"

samtools index $out/$d/$d.sorted.markdu.bam \\
&& echo "bam index done!"


rm $out/$d/$d.sorted.bam

END_BAM
    print $ou2 "$bam\n";
    print $sou "subbam$d=\$(sbatch --dependency=afterok:\$jidm$d -c 4 --mem=10G $d.sort_markdu_bam.sh)\njidb$d=\$(echo \$subbam$d | awk -F' ' '{print \$4}')\n";


    for my $chr (@Chrs) {
            open my $ou3, ">", "$shdir/$d.$chr.hapcall.sh" or die;
            my $jobchr = $d.$chr;
        my $hapcall = <<"END_CALL";
#!/bin/bash
#SBATCH -J $jobchr
#SBATCH -o $jobchr.out
#SBATCH -e $jobchr.err
module purge

$gatk

gatk --java-options "-Xmx8g" HaplotypeCaller \\
-R $genome \\
-L $chr \\
-ploidy 2 \\
--emit-ref-confidence GVCF \\
-I $out/$d/$d.sorted.markdu.bam \\
-O $out/$d/$d.$chr.g.vcf.gz \\
&& echo "call $d $chr gvcf done!"

gatk --java-options "-Xmx8g" IndexFeatureFile \\
-I $out/$d/$d.$chr.g.vcf.gz \\
&& echo "index $d $chr gvcf done!"

END_CALL
    print $ou3 "$hapcall\n";
    print $sou "chr=\$(sbatch --dependency=afterok:\$jidb$d -c 4 --mem=10G $d.$chr.hapcall.sh)\n";
    }
}

