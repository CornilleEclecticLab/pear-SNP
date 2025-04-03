#!/usr/bin/env bash

#SBATCH -J s02.match_genes_gff_with_annotation.sh
#SBATCH -o s02.match_genes_gff_with_annotation.sh.%J.out
#SBATCH -e s02.match_genes_gff_with_annotation.sh.%J.err

# @File     :   s02.match_genes_gff_with_annotation.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/04/02 16:23:33
#Description:
    # 

source load s00.load_bedtools.sh

for bed in $(ls ../output/*/NLGenomeSweeper/Final_candidates.bed); do
    grep "NL" $bed > $bed.NLs.bed

    genome_name=$(basename $(dirname $(dirname $bed)))

    gff="../input/${genome_name}.chr_list.gene.gff"

    bedtools intersect -a $gff -b $bed.NLs.bed -wa -wb -F 1.0 > $bed.NLs.gff

    cat $bed.NLs.gff \
    | awk -F'[;=\t ]' '{for(i=1;i<=NF;i++) if($i=="Accession" || $i=="Name") print $(i+1) "\t" $NF}' \
    | sort -u \
    > $bed.NLs.genes.txt

done

