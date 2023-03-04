#!/bin/bash
# Adaptation by Yuqi NIE on 2023/01/26

#SBATCH -J s04.filter_syno
#SBATCH -o s04.filter_syno.out
#SBATCH -e s04.filter_syno.err

module purge
module load bioinfo/samtools-1.14

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

perl filter_by_anno.pl \
$OUTPUT/pear.Combine_Chr.geno20_maf005.anno.vcf.gz \
$OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.vcf

bgzip $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.vcf && \
tabix -p vcf $OUTPUT/pear.Combine_Chr.geno20_maf005.anno.syno.vcf.gz
