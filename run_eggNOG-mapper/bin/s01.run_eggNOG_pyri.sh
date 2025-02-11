#!/bin/bash

#SBATCH -J eggNOG.sh
#SBATCH -o s01.run_eggNOG_pyri.sh.%J.out
#SBATCH -e s01.run_eggNOG_pyri.sh.%J.err
#SBATCH -c 12
#SBATCH --mem=60G
#SBATCH -p long

# @File     :   s01.run_eggNOG.sh
# @Version  :   1.2.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/12/05 15:50:35
#Description:
    # 

module load eggnog-mapper/2.1.12

TMP=../tmp_eggnog-mapper
mkdir -p $TMP

# run eggnog-mapper
emapper.py \
    -m diamond \
    -i ../input/GWHBAOS00000000.Protein.faa \
    --itype proteins \
    -o ../output/GWHBAOS00000000 \
    --dbmem \
    --dmnd_db ../database/Viridiplantae.dmnd \
    --temp_dir $TMP \
    --data_dir ../database \
    --override \
    --cpu 0
