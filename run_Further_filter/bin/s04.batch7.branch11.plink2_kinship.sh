#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship.out
#SBATCH -e s04.kinship.err
module purge
module load bioinfo/plink2-v2.0_alpha2
module load bioinfo/plink-v1.90b5.3

mkdir -p /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/s04.clone_filtering

plink --vcf /work/zruilin/pear/run_Population_genetic_filter/output/branch7.batch11.whole_pear_keep_all_pear/whole_pear.keep_all_samples.Combine_Chr.geno20_maf005.vcf.gz \
--make-bed \
--const-fid \
--out /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/whole_pear.keep_all_samples.Combine_Chr.geno20_maf005
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr




# Do not change --bfile to --vcf 
	# If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
	# Invalid chromosome code 'Pp01' on line 240 of --vcf file.
	# (Use --allow-extra-chr to force it to be accepted.)
plink2 --bfile /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/whole_pear.keep_all_samples.Combine_Chr.geno20_maf005 \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/s04.clone_filtering/without_clone0.354

plink2 --bfile /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/whole_pear.keep_all_samples.Combine_Chr.geno20_maf005 \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.125 \
--make-bed \
--out /work/zruilin/pear/run_Further_filter/output/branch7.batch11.whole_pear_keep_all_pear/s04.clone_filtering/without_2nd-degree_relation0.125


