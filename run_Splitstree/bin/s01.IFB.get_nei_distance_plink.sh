#!/bin/bash
#SBATCH -J s01.get_nei_plink
#SBATCH -o s01.get_nei_plink.%J.out
#SBATCH -e s01.get_nei_plink.%J.err
module purge
module load plink/1.90b6.18


WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_Splitstree/bin"
INPUT="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/branch8.batch15.pear_Jul2023_noclone/pear_Jul2023_noclone.Combine_Chr.geno20_maf005.anno.syno.thin8k"

plink -bfile $INPUT --allow-extra-chr --distance 1-ibs square
