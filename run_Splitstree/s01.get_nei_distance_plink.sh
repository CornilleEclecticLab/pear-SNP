#!/bin/bash
#SBATCH -J s01.get_nei_plink
#SBATCH -o s01.get_nei_plink.out
#SBATCH -e s01.get_nei_plink.err
module purge
module load bioinfo/plink-v1.90b5.3


WORKDIR="/home/zruilin/work/Peach/run_Splitstree"
INPUT="$WORKDIR/input/Combine_Chr.geno20_maf005.anno.syno.thin8k"

plink -bfile $INPUT --allow-extra-chr --distance 1-ibs square
