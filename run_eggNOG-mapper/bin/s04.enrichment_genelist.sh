#!/bin/bash

#SBATCH -J s04.enrichment_genelist.sh
#SBATCH -o s04.enrichment_genelist.sh.%J.out
#SBATCH -e s04.enrichment_genelist.sh.%J.err


# Rscipt s04.x1.enrichment_and_plot.R <DB_NAME[comm | pyri]> <gene_list_file> <output_dir>

mkdir -p ../output/s04.pyri
for i in ../../run_OmegaPlus/output/s08.summary_table/specific.positive_selection_genes_list.*.txt \
        ../../run_OmegaPlus/output/s08.summary_table/interest_common_pop.positive_selection_genes_list.*.txt
    do Rscript s04.x1.enrichment_and_plot.R \
        pyri \
        $i \
        ../output/s04.pyri
done


mkdir -p ../output/s04.comm
for j in ../../run_OmegaPlus_refer_on_communis/output/s08.summary_table/interest_common_pop.positive_selection_genes_list.*.txt \
         ../../run_OmegaPlus_refer_on_communis/output/s08.summary_table/specific.positive_selection_genes_list.*.txt
    do Rscript s04.x1.enrichment_and_plot.R \
        comm \
        $j \
        ../output/s04.comm/
 done