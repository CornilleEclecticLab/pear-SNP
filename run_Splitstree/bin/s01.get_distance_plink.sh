#!/bin/bash
#SBATCH -J s01.get_distance_plink
#SBATCH -o s01.get_distance_plink.out
#SBATCH -e s01.get_distance_plink.err
module purge
module load bioinfo/plink-v1.90b5.3


WORKDIR="/work/zruilin/pear/run_Splitstree"
INPUT="/work/zruilin/pear/run_Population_genetic_filter/output/batch1.reduce3EupropeWildPear/pear.Combine_Chr.geno20_maf005.anno.syno.thin8k"

plink -bfile $INPUT --allow-extra-chr --distance 1-ibs square
