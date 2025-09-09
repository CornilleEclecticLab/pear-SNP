#!/bin/bash
# By Xilong CHEN
# Create date: 2025-06-10
# Contact: chen_xilong@outlook.com
# Modified by Yuqi on 2025-09-01

#SBATCH -c 40
#SBATCH -e astral4.5pop.%j.err
#SBATCH -o astral4.5pop.%j.out
#SBATCH --mem=500G


module purge
# module load parallel/20200822
module load parallel/20190322

# source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
# conda activate Aphid

module load conda
. activate aster

mkdir -p ../output/ASTER.5pop/astral4.5pop.RES

astral4 \
    -t 40 \
    -i ../output/ASTER.5pop/astral4.5pop.RES/auto.5pop.trees \
    -o ../output/ASTER.5pop/astral4.5pop.RES/auto5pop.sptree.tre \
    -a ../output/ASTER.5pop/astral4.5pop.RES/ussuWest_AB_group.txt \
    --root ussu

