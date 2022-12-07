#!/usr/bin/env bash


# @File     :   s04.plink2.kinship.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2022/12/06 16:08:24
#Description:
	# 
module purge
module load bioinfo/plink2-v2.0_alpha2
module load bioinfo/plink-v1.90b5.3



# Transform vcf to bfile
plink --vcf ../output/Pear.Zhang_2021.test_rand.vcf.gz \
--make-bed \
--const-fid \
--out ../output/Pear.Zhang2021.test_rand
--set-missing-var-ids @:# \
--keep-allele-order \
--allow-extra-chr


plink2 --bfile ../output/Pear.Zhang_2021.test_rand \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.354 \
--make-bed \
--out ../output/without_clone/Pear.cutoff0.354


plink2 --bfile ../output/Pear.Zhang_2021.test_rand \
--allow-extra-chr \
--make-king-table \
--king-cutoff 0.400 \
--make-bed \
--out ../output/without_clone/Pear.cutoff0.400