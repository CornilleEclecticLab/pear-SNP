#!/bin/bash
#SBATCH -J s01.run_Dtest
#SBATCH -o s01.run_Dtest.%J.out
#SBATCH -e s01.run_Dtest.%J.err
#SBATCH -c 2
#SBATCH --mem=4G

source s00.load_environment.sh

Dsuite  Dtrios \
    -o $WORK_DIR/output/output_full_tree_full_list \
    -t $WORK_DIR/input/s01.full_tree.nwk  \
    $WORK_DIR/input/s01.input.vcf  \
    $WORK_DIR/input/s01.id_species.txt 

Dsuite  Dtrios \
    -o $WORK_DIR/output/output_full_tree_Western_list \
    -t $WORK_DIR/input/s01.full_tree.nwk  \
    $WORK_DIR/input/s01.input.vcf  \
    $WORK_DIR/input/s01.id_species.Western.txt 

Dsuite  Dtrios \
    -o $WORK_DIR/output/output_full_tree_Eastern_list \
    -t $WORK_DIR/input/s01.full_tree.nwk  \
    $WORK_DIR/input/s01.input.vcf  \
    $WORK_DIR/input/s01.id_species.Eastern.txt 

Dsuite  Dtrios \
    -o $WORK_DIR/output/output_Eastern_tree_Eastern_list \
    -t $WORK_DIR/input/s01.Eastern_tree.nwk  \
    $WORK_DIR/input/s01.input.vcf  \
    $WORK_DIR/input/s01.id_species.Eastern.txt 

Dsuite  Dtrios \
    -o $WORK_DIR/output/output_Western_tree_Western_list \
    -t $WORK_DIR/input/s01.Western_tree.nwk  \
    $WORK_DIR/input/s01.input.vcf  \
    $WORK_DIR/input/s01.id_species.Western.txt
