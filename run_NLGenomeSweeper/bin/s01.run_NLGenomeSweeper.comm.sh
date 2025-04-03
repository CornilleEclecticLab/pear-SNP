#!/bin/bash

#SBATCH -J NLGcomm
#SBATCH -o s01.run_NLGenomeSweeper.comm.sh.%J.out
#SBATCH -e s01.run_NLGenomeSweeper.comm.sh.%J.err
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -p long

# @File     :   s01.run_NLGenomeSweeper.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/03/20 14:50:59
#Description:
    # 

module load conda
source activate NLGenomeSweeper

outdir=../output/PyrusCommunis_BartlettDHv2.0

mkdir -p $outdir

NLGenomeSweeper \
    -genome ../input/PyrusCommunis_BartlettDHv2.0.fasta \
    -t 8 \
    -outdir $outdir
