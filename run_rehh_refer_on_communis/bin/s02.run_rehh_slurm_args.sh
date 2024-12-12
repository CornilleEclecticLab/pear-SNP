#!/usr/bin/env bash

#SBATCH -J rehh_args
#SBATCH -o s02.run_rehh_slurm_args.sh.%J.out
#SBATCH -e s02.run_rehh_slurm_args.sh.%J.err
#SBATCH -c 18
#SBATCH --mem=64G

# @File     :   s02.run_rehh_slurm.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/12/10 22:06:59
#Description:
    # 

# load R and parameters
source s00.config.sh

# CPU=$(($SLURM_CPUS_PER_TASK - 1))

Rscript s02.x1.rehh.args.R \
    -s $SAMPLE \
    -o $OUTGROUP \
    -v $VCF_FOLDER \

