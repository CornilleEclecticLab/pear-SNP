#!/bin/bash
#SBATCH -J s02.run_Dtest_by_chr.root_P.ussu_Western_tree_Western_list.geno20
#SBATCH -o s02.run_Dtest_by_chr.root_P.ussu_Western_tree_Western_list.geno20.%J.out
#SBATCH -e s02.run_Dtest_by_chr.root_P.ussu_Western_tree_Western_list.geno20.%J.err
#SBATCH -c 20
#SBATCH --mem=200G


source s00.load_environment.sh
RUN_NAME="root_P.ussu_Western_tree_Western_list.geno20"

# DtriosParallel \
#     --run-name "${RUN_NAME}" \
#     --keep-intermediate \
#     --tree "$WORK_DIR/output/${RUN_NAME}/root_ussu_Western.nwk" \
#     "$WORK_DIR/output/${RUN_NAME}/s02.id_species.Western.outgroup_cauc.txt" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000076.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000085.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000099.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000128.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000158.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000163.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000172.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000224.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000274.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000335.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000352.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000356.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000365.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000381.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000386.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000424.choosing.vcf.gz" \
#     "$WORK_DIR/output/s01.split_vcf/GWHBAOS00000425.choosing.vcf.gz" 



source s00.load_ruby.sh

for mode in tree Dmin BBAA; do
    ruby plot_d.NOGIT.rb \
        "$WORK_DIR/output/$RUN_NAME/DTparallel_s02.id_species.Western.outgroup_cauc_${RUN_NAME}_combined_${mode}.txt" \
        "$WORK_DIR/output/$RUN_NAME/population_order.txt" \
        0.7 \
        $WORK_DIR/output/$RUN_NAME/plot.${mode}_D.cutoff07.svg

    ruby plot_f4ratio.NOGIT.rb \
        "$WORK_DIR/output/$RUN_NAME/DTparallel_s02.id_species.Western.outgroup_cauc_${RUN_NAME}_combined_${mode}.txt" \
        "$WORK_DIR/output/$RUN_NAME/population_order.txt" \
        0.2 \
        $WORK_DIR/output/$RUN_NAME/plot.${mode}_f4ratio.cutoff02.svg
done
