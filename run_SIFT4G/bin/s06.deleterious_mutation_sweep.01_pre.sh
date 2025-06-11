#!/bin/bash


#SBATCH -J s06.deleterious_mutation_sweep.01_pre.sh
#SBATCH -o s06.deleterious_mutation_sweep.01_pre.sh.%J.out
#SBATCH -e s06.deleterious_mutation_sweep.01_pre.sh.%J.err
#SBATCH -c 1
#SBATCH --mem=48G


module load conda
source activate intervaltree

## positive selection
python3 s06.x1.dm_extract_sweep.py \
    ../output/s06.dm_sweep \
    selection \
    ../input/selection_bed/selection.outliers.betu.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.cauc.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.comm_Dessert.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.comm_Perry.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.pyra.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.pyri_JP.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.Sand_CN-SE.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.Sand_CN-SW.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.ussu.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.White.filtered.sorted.merged.bed

## control
python3 s06.x1.dm_extract_sweep.py \
    ../output/s06.dm_sweep \
    control \
    ../input/selection_bed/selection.outliers.betu.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.cauc.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.comm_Dessert.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.comm_Perry.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.pyra.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.pyri_JP.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.Sand_CN-SE.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.Sand_CN-SW.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.ussu.filtered.sorted.merged.bed \
    ../input/selection_bed/selection.outliers.White.filtered.sorted.merged.bed

## all_masked
python3 s06.x1.dm_extract_sweep.py \
    ../output/s06.dm_sweep \
    all_masked 


## Combine the output files into one file
head -n1 ../output/s06.dm_sweep/control.pop_mutation_burden.txt > ../output/s06.dm_sweep/region.pop_mutation_burden.Rdata.txt

for file in ../output/s06.dm_sweep/*.pop_mutation_burden.txt; do 
    if [[ "$file" != "../output/s06.dm_sweep/region.pop_mutation_burden.Rdata.txt" ]]; then
        tail -n +2 "$file" >> ../output/s06.dm_sweep/region.pop_mutation_burden.Rdata.txt
    fi
done


echo End Time :
date
