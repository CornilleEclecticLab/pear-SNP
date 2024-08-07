#!/usr/bin/env bash

#SBATCH -J callable_sites_mask.sh
#SBATCH -o callable_sites_mask.sh.%J.out
#SBATCH -e callable_sites_mask.sh.%J.err
#SBATCH -c 4
#SBATCH --mem=360G


# @File     :   callable_sites_mask.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/28 18:33:19
#Description:
    #

module load conda 
source activate bedops
module load bedtools/2.30.0

cd ../input/allsites_vcf_files

for i in *.vcf.gz ; do 
    echo "Reading" $i
    zcat $i | grep -v "#" | awk '{print $1"\t"$2"\t"$2}' | bedops --merge - > ${i%.vcf.gz}.callable_sites.bed &
    zcat $i | vcf2bed > ${i%.vcf.gz}.mask_sites.bed &
done

wait

cat *.mask_sites.bed >  all.mask_sites.bed

cat *.callable_sites.bed > all.callable_sites.bed

bedtools merge -i all.mask_sites.bed > all.mask_sites.merged.bed
