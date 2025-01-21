#!/usr/bin/env bash

#SBATCH -J s03.run_bedtools_intersect.sh
#SBATCH -o s03.run_bedtools_intersect.sh.%J.out
#SBATCH -e s03.run_bedtools_intersect.sh.%J.err

# @File     :   s03.run_bedtools_intersect.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/01/20 17:12:23
#Description:
    # 

source s00.load_bedtools.sh

# sort the bed files
for bed in $(ls ../output/*.gene_density_pass.bed); do
    base=$(basename $bed)
    prefix=$(echo $base | sed 's/s02.R.//'| sed 's/.gene_density_pass.bed//')

    awk '{print $1"\t"$2"\t"$3}' ${bed} | bedtools sort > ${bed}.sorted.bed
    # ls ${bed}.sorted.bed
    # ls ../input/${prefix}.genome.genmap_pass.chr_list.sorted.bed
    bedtools intersect \
        -a ${bed}.sorted.bed \
        -b ../input/${prefix}.genome.genmap_pass.chr_list.sorted.bed \
    | bedtools sort \
    | bedtools merge \
    > ${bed}.genmap_density_pass_intersect.bed

    while read -r chr;do
        grep ${chr} ${bed}.genmap_density_pass_intersect.bed \
        > ../output/${chr}.genmap_density_pass_intersect.bed
    done < ../input/${prefix}.chr_list.txt
done
