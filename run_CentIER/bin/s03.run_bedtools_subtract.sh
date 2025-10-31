#!/usr/bin/env bash

#SBATCH -J s03.run_bedtools_subtract.sh
#SBATCH -o s03.run_bedtools_subtract.sh.%J.out
#SBATCH -e s03.run_bedtools_subtract.sh.%J.err

# @File     :   s03.run_bedtools_subtract.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/01/29 13:59:18
#Description:
    # 

source s00.load_bedtools.sh

# sort the bed files
suffix=".chr_list.fasta_centromere_range.txt"

for bed in $(ls ../output/s02.generate_run_CentIER/*/*${suffix}); do
    base=$(basename $bed)
    prefix=$(echo $base | sed "s/$suffix//")

    awk '{print $1"\t"$2-1"\t"$3}' ${bed} \
    | awk -v OFS='\t'  'NR==FNR {map[$2] = $1; next} $1 in map { $1 = map[$1];print } 1' ../input/${prefix}.chrID_map.txt ${bed}.bed - \
    | bedootls sort \
    > ${bed}.sorted.bed

    bedtools subtract \
        -a ../input/${prefix}.genome.genmap_pass.chr_list.sorted.bed \
        -b ${bed}.sorted.bed \
    | bedtools sort \
    | bedtools merge \
    > ${bed}.genmap_pass_subtract_centier.bed

done
