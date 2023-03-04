#!/bin/bash
# By Yuqi NIE on 2022/11/03

#SBATCH -J s05a.LD_stat
#SBATCH -o s05a.LD_stat.out
#SBATCH -e s05a.LD_stat.err

module purge
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15

WORKDIR="/work/zruilin/Peach/run_Further_filter"
INPUT="$WORKDIR/input"
OUTPUT="$WORKDIR/output"

vcftools --vcf $INPUT/vcfTest.test.vcf \
--hap-r2 \
--ld-window-bp 500 \
--out $OUTPUT/s05a.ld_windows_500vcftools \

