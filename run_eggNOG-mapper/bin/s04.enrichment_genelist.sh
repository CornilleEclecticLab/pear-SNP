#!/bin/bash

#SBATCH -J s04.enrichment_genelist.sh
#SBATCH -o s04.enrichment_genelist.sh.%J.out
#SBATCH -e s04.enrichment_genelist.sh.%J.err
#SBATCH -a 1-2

module load r/4.3.1


# Rscipt s04.x1.enrichment_and_plot.R <DB_NAME[comm | pyri]> <gene_list_file> <output_dir>

mkdir -p ../output/s04.pyri_and_comm

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
for i in ../../run_OmegaPlus/output/s08.a1.summary_gene_list/specific.positive_selection_genes_list.*.txt \
         ../../run_OmegaPlus/output/s08.a1.summary_gene_list/shared_any_cultivar.positive_selection_gene_list.cultivar.txt \
         ../../run_OmegaPlus/output/s08.a1.summary_gene_list/shared_per_wild_cultivar.positive_selection_gene_list.pash.txt
    do Rscript s04.x1.enrichment_and_plot.R \
        pyri \
        $i \
        ../output/s04.pyri_and_comm/
    done
fi


if [ $SLURM_ARRAY_TASK_ID -eq 2 ]; then
for j in ../../run_OmegaPlus_refer_on_communis/output/s08.a1.summary_gene_list/specific.positive_selection_genes_list.*.txt \
         ../../run_OmegaPlus_refer_on_communis/output/s08.a1.summary_gene_list/shared_any_cultivar.positive_selection_gene_list.cultivar.txt \
         ../../run_OmegaPlus_refer_on_communis/output/s08.a1.summary_gene_list/shared_any_wild_cultivar.positive_selection_gene_list.wild_cultivar.txt
    do Rscript s04.x1.enrichment_and_plot.R \
        comm \
        $j \
        ../output/s04.pyri_and_comm/
    done
fi