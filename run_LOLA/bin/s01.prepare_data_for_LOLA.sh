#!/bin/bash

#SBATCH -J s01.prepare_data_for_LOLA.sh
#SBATCH -o s01.prepare_data_for_LOLA.sh.%J.out
#SBATCH -e s01.prepare_data_for_LOLA.sh.%J.err
#SBATCH --mem=50G
#SBATCH -c 4


# @File     :   s01.prepare_data_for_LOLA.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/08/08 17:59:09
#Description:
    #

module purge
. s00.load_modules.sh

set -euo pipefail

path=$(dirname $(pwd))
mkdir -p $path/input/bed
mkdir -p $path/input/gff



# The two scripts run one time is OK for both PPY and PCOM

# Selective gene ->  GFF format is OK
python3 s01.x2.get_selection_gene_gff.py

# Get masked genes regions -> bed format
python3 s01.x1.mask_gene.py


samples=(PPY PCOM)
for sample in "${samples[@]}"; do

cd $path/bin

. s00.load_samtools.sh
samtools faidx $path/input/genomes/${sample}.fasta
. s00.unload_samtools.sh # conda env of samtools would cause bash error libtinfo.so.6

# Get control genes regions: masked genes - selection genes -> bed format
bedtools intersect -v \
    -a $path/input/bed/${sample}.gene_masked.bed \
    -b $path/input/gff/${sample}.selection.gff \
    > $path/input/bed/${sample}.control.bed


bedtools intersect -v \
    -a $path/input/gff/${sample}.gene_masked.gff \
    -b $path/input/gff/${sample}.selection.gff \
    > $path/input/gff/${sample}.control.gff

cd $path/input/bed

    for slop in 2 10; do
        # slop
        ## If use gff, can get upstream 2kb, rather -b both stream
        bedtools slop -i ${sample}.gene.bed -g $path/input/genomes/${sample}.fasta.fai -b ${slop}000 > ${sample}.gene_${slop}kb.bed
        bedtools slop -i ../gff/${sample}.selection.gff -g ../genomes/${sample}.fasta.fai -b ${slop}000 > ../gff/${sample}.selection_${slop}kb.gff
        bedtools slop -i ${sample}.control.bed -g ../genomes/${sample}.fasta.fai -b ${slop}000 > ${sample}.control_${slop}kb.bed

        # sort
        bedtools sort -i ${sample}.gene_${slop}kb.bed > ${sample}.gene_${slop}kb.sorted.bed
        bedtools sort -i ../gff/${sample}.selection_${slop}kb.gff > ../gff/${sample}.selection_${slop}kb.sorted.gff
        bedtools sort -i ${sample}.control_${slop}kb.bed > ${sample}.control_${slop}kb.sorted.bed


        # merge, gff will be bed
        bedtools merge -i ${sample}.gene_${slop}kb.sorted.bed > ${sample}.gene_${slop}kb.sorted.merged.bed
        bedtools merge -i ../gff/${sample}.selection_${slop}kb.sorted.gff > ../bed/${sample}.selection_${slop}kb.sorted.merged.bed
        bedtools merge -i ${sample}.control_${slop}kb.sorted.bed > ${sample}.control_${slop}kb.sorted.merged.bed

        # rm temp files
        # rm ${sample}.gene_${slop}kb.bed
        # rm ../gff/${sample}.selection_${slop}kb.gff
        # rm ${sample}.control_${slop}kb.bed

        # rm  ${sample}.gene_${slop}kb.sorted.bed
        # rm ../gff/${sample}.selection_${slop}kb.sorted.gff
        # rm ${sample}.control_${slop}kb.sorted.bed
    done

cd $path/input/TE_annotation_URGI/${sample}

if [ $sample == 'PPY' ]; then
    sample_number='PPY3'  # PPY is PPY3
else
    sample_number=$sample
fi

echo "Building LOLA DB... " $(pwd)
cat ${sample_number}_TEannotGr2_GFF3chr/*.gff3 > ${sample}_TE_annotation.gff3
grep '\smatch' ${sample}_TE_annotation.gff3 | grep -v match_part > ${sample}_TE_annotation.filtered.gff3

if [ $sample == "PPY" ]; then
    python3 $path/bin/s03.x2.pure_chrom_in_TE_gff.py ${sample}_TE_annotation.filtered.gff3
fi

{
    perl $path/bin/s01.x3.REPEAT2LOLA.pl -gff ${sample}_TE_annotation.filtered.gff3
    rm -rf TE_DB/${sample}
    mkdir -p TE_DB/${sample}
    mv regions TE_DB/${sample}/
    echo "LOLA DB for ${sample} completed"
} &

echo "Transforming GFF to TSV..."
# Transfer GFF to TSV for s03.tegrip.R
python3 $path/bin/s03.x1.transform_gff_to_tsv.py \
    --gff $path/input/gff/${sample}.control.gff \
    --mode gene \
    --output_tsv $path/input/${sample}.control.gene_annotation.tsv \
    --classif control

python3 $path/bin/s03.x1.transform_gff_to_tsv.py \
    --gff $path/input/gff/${sample}.gene_masked.gff \
    --mode gene \
    --output_tsv $path/input/${sample}.gene_masked.gene_annotation.tsv \
    --classif gene_masked

python3 $path/bin/s03.x1.transform_gff_to_tsv.py \
    --gff $path/input/gff/${sample}.selection.gff \
    --mode gene \
    --output_tsv $path/input/${sample}.selection.gene_annotation.tsv \
    --classif selection

echo "Transforming TE GFF to TSV"
python3 $path/bin/s03.x1.transform_gff_to_tsv.py \
    --gff ${sample}_TE_annotation.filtered.gff3 \
    --mode TE \
    --output_tsv $path/input/${sample}.TE_annotation.tsv \
    --classif ${sample_number}_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif &
done

wait
