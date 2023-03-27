#!/bin/bash
# Adaptation by Yuqi NIE on 2023/01/26

#SBATCH -J s04.filter_syno
#SBATCH -o s04.filter_syno.out
#SBATCH -e s04.filter_syno.err

module purge
module load bcftools/1.14

WORKDIR="/shared/ifbstor1/projects/pear_snp2/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.whole_pear"
PREFIX="whole_pear"

perl filter_by_anno.pl \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf.gz \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf

bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf.gz
