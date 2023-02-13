#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship.out
#SBATCH -e s04.kinship.err
module purge
module load bioinfo/plink2-v2.0_alpha2
module load bioinfo/plink-v1.90b5.3

plink --vcf /home/zruilin/work/Peach/run_Further_filter/output/vcfTest.test.rand500.vcf \
--make-bed \
--const-fid \
--out /home/zruilin/work/Peach/run_Further_filter/output/vcfTest.test.rand500 \
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr


# Do not change --bfile to --vcf 
	# If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
	# Invalid chromosome code 'Pp01' on line 240 of --vcf file.
	# (Use --allow-extra-chr to force it to be accepted.)
plink2 --bfile /home/zruilin/work/Peach/run_Further_filter/output/vcfTest.test.rand500 \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out /home/zruilin/work/Peach/run_Further_filter/output/s04.clone_filtering/without_clone
