#!/bin/bash

#SBATCH -J s03.snpEff_annotaion
#SBATCH -o s03.snpEff_annotaion.%j.out
#SBATCH -e s03.snpEff_annotaion.%j.err

module purge
module load bcftools/1.14
module load samtools/1.14
module load java-jdk/11.0.9.1


WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch8.batch15.pear_Jul2023_noclone"
PREFIX="pear_Jul2023_noclone"

tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz

java -Xmx8g -jar $WORKDIR/bin/snpEff/snpEff.jar Pyr.cuiguan \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz > \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf

mv snpEff_* $OUTPUT/

bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf.gz
