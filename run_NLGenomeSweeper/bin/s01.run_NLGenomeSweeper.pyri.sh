#!/bin/bash

#SBATCH -J NLGpyri
#SBATCH -o s01.run_NLGenomeSweeper.pyri.sh.%J.out
#SBATCH -e s01.run_NLGenomeSweeper.pyri.sh.%J.err
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

outdir=../output/GWHBAOS00000000

mkdir -p $outdir

NLGenomeSweeper \
    -genome ../input/GWHBAOS00000000.fasta \
    -t 8 \
    -outdir $outdir
