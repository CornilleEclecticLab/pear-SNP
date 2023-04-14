#!/bin/bash

#SBATCH -J s03.snpEff_annotaion
#SBATCH -o s03.snpEff_annotaion.out
#SBATCH -e s03.snpEff_annotaion.err

module purge
module load bioinfo/bcftools-1.14
module load bioinfo/samtools-1.14
module load bioinfo/Java15.0.1

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.batch12.kinship_base_on_branch7_batch11_remove_clone_low_quality"
PREFIX="whole_pear"

tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz

java -Xmx8g -jar $WORKDIR/bin/snpEff/snpEff.jar Pyr.cuiguan \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.vcf.gz > \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf

mv snpEff_* $OUTPUT/

bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf.gz
