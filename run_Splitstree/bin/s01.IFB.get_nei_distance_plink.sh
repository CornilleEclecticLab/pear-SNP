#!/bin/bash
#SBATCH -J s01.get_nei_plink
#SBATCH -o s01.get_nei_plink.%J.out
#SBATCH -e s01.get_nei_plink.%J.err
module purge
module load plink/1.90b6.18

source s00.config.sh

plink -bfile "$INPUT" --allow-extra-chr --distance 1-ibs square

python3 s01.x1.format_mdist.py
