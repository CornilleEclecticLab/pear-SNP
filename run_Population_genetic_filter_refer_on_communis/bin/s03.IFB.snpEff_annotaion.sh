#!/bin/bash

#SBATCH -J s03.snpEff_annotaion
#SBATCH -o s03.snpEff_annotaion.NOGIT.%j.out
#SBATCH -e s03.snpEff_annotaion.NOGIT.%j.err

module purge
module load bcftools/1.14
module load samtools/1.14
module load java-jdk/11.0.9.1
source s00.config.sh

tabix -p vcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.vcf.gz

java -Xmx8g -jar $WORKDIR/bin/snpEff.NOGIT./snpEff.jar  PyrusCommunis_BartlettDHv2.0 \
$OUTPUT/$PREFIX.Combine_chr.geno20_maf005.vcf.gz \
    > $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.vcf

mv snpEff_* $OUTPUT/

bgzip $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.vcf.gz
