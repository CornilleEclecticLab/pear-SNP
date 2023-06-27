#!/bin/bash

#SBATCH -J kinship
#SBATCH -o s04.kinship.%j.out
#SBATCH -e s04.kinship.%j.err
module purge
module load bioinfo/plink2-v2.0_alpha2
module load bioinfo/plink-v1.90b5.3


OUTDIR=/work/zruilin/pear/run_Further_filter/output/branch6.loquat_and_loquat2/s04.clone_filtering
BED="/work/zruilin/pear/run_Further_filter/output/branch6.loquat_and_loquat2/loquat.Combine_Chr"

mkdir -p $OUTDIR

plink --vcf /work/zruilin/pear/run_Population_genetic_filter/output/branch6.loquat_and_loquat2/loquat.Combine_Chr.vcf.gz  \
--make-bed \
--const-fid \
--out $BED \
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr




# Do not change --bfile to --vcf 
	# If you only have vcf files, use the upwards to transform them into bfile format. 
#Add --allow-extra-chr 
	# Invalid chromosome code 'Pp01' on line 240 of --vcf file.
	# (Use --allow-extra-chr to force it to be accepted.)
plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out $OUTDIR/without_clone0.354

plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.125 \
--make-bed \
--out $OUTDIR/without_2nd-degree_relation0.125

plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.177 \
--make-bed \
--out $OUTDIR/without_2nd-degree_relation0.177

plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.445 \
--make-bed \
--out $OUTDIR/high_cut_off0.445

plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.446 \
--make-bed \
--out $OUTDIR/high_cut_off0.446

plink2 --bfile $BED \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.447 \
--make-bed \
--out $OUTDIR/high_cut_off0.447
