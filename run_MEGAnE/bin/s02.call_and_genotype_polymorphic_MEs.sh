#!/usr/bin/env bash

#SBATCH -J s02.call_and_genotype_polymorphic_MEs.sh
#SBATCH -o s02.call_and_genotype_polymorphic_MEs.sh.%J.out
#SBATCH -e s02.call_and_genotype_polymorphic_MEs.sh.%J.err

# @File     :   s02.call_and_genotype_polymorphic_MEs.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 18:03:22
#Description:
    # 

source s00.config.sh


# In the case of CRAM file mapping to GRCh38-related genome (e.g. GRCh38DH, hg38)
singularity exec ${sif} call_genotype_Cuiguan \
-i /path/to/input.cram \
-fa genome_file.fasta \
-mk ./megane_kmer_set/reference_Pyrifolia_Cuiguan.mk \
-outdir MEGAnE_result_test \
-sample_name test_sample \
-p 4