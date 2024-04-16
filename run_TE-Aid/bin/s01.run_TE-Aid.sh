#!/usr/bin/env bash

#SBATCH -J s01.run_TE-Aid.sh
#SBATCH -o s01.run_TE-Aid.sh.%J.out
#SBATCH -e s01.run_TE-Aid.sh.%J.err

# @File     :   s01.run_TE-Aid.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/08 22:51:12
#Description:
    # 


module load conda
source activate TE_AID

# Run TE-Aid
./TE-Aid/TE-Aid \
    --query ../input/reduced_Inpactor2_library_c10.fasta \
    --genome ../input/GWHBAOS00000000.genome.fasta \
    --output ../output \
    --emboss-dotmatcher
