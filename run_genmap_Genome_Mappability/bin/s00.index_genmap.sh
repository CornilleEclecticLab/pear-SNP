#!/usr/bin/env bash

# @File     :   s00.index_genmap.bash
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/07/19 17:49:29
#Description:
	# 


#SBATCH -J index_genmap
#SBATCH -e s00.index_genmap.NOGIT.%J.err
#SBATCH -o s00.index_genmap.NOGIT.%J.out
#SBATCH -c 1
#SBATCH --mem=24G


module load conda

source activate ~/work/conda/env/genmap

# genmap index -v -F ../input/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.fasta -I ../input/index_genmap
genmap index -v -F ../input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta -I ../input/index_genmap_Pyrus_communis
