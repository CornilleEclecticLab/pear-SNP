#!/bin/bash
# Adaptation by Yuqi NIE on 2023/01/26

#SBATCH -J s04.filter_syno
#SBATCH -o s04.filter_syno.%j.out
#SBATCH -e s04.filter_syno.%j.err

module purge
module load bioinfo/samtools-1.14

WORKDIR="/work/zruilin/pear/run_Population_genetic_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output/branch7.batch12.kinship_base_on_branch7_batch11_remove_clone_low_quality"
PREFIX="whole_pear"

perl filter_by_anno.pl \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.vcf.gz \
$OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf

bgzip $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf && \
tabix -p vcf $OUTPUT/$PREFIX.Combine_Chr.geno20_maf005.anno.syno.vcf.gz
