#!/bin/bash
# Adaptation by Yuqi NIE on 2023/01/26

#SBATCH -J s04.filter_syno
#SBATCH -o s04.filter_syno.NOGIT.%j.out
#SBATCH -e s04.filter_syno.NOGIT.%j.err

module purge
module load bcftools/1.14
source s00.config.sh

perl filter_by_anno.pl \
$OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.vcf.gz \
$OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf

bgzip $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_chr.geno20_maf005.anno.syno.vcf.gz
