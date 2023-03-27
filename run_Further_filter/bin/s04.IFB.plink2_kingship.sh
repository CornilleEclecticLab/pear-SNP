#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship.%J.out
#SBATCH -e s04.kinship.%J.err

module purge

module load plink2
# PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)
module load plink
# PLINK v1.90b5 64-bit (14 Nov 2017)

plink --vcf /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/whole_pear.combine.test_rand.vcf.gz \
--make-bed \
--const-fid \
--out /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/whole_pear.combine.test.rand \
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr


mkdir -p /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/s04.clone_filtering

# Do not change --bfile to --vcf 
	# If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
	# Invalid chromosome code 'Pp01' on line 240 of --vcf file.
	# (Use --allow-extra-chr to force it to be accepted.)
plink2 --bfile /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/whole_pear.combine.test.rand \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/s04.clone_filtering/without_clone0.354

plink2 --bfile /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/whole_pear.combine.test.rand \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.125 \
--make-bed \
--out /shared/ifbstor1/projects/pear_snp2/pear/run_Further_filter/output/whole_pear/s04.clone_filtering/without_2nd-degree_relation0.125


