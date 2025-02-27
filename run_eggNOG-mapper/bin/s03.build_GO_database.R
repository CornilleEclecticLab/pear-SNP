# By Xilong CHEN
# Create date: 2023-01-20
# Contact: chen_xilong@outlook.com

# Modified by Yuqi NIE on 2025-02-11

# Install AnnotationForge
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationForge")

# Set Working Directory to source file folder
rm(list = ls())
library(AnnotationForge)
library(tibble)
library(jsonlite)
library(stringr)

output_dir <- "../output/s03.build_GO_KEGG_database/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Genome 1
gene_info <- read.table("../output/s02.prepare_GO_KEGG/PyrusCommunis_BartlettDHv2.0.GO_GID.table.txt", head =T)
gene2go <- read.table("../output/s02.prepare_GO_KEGG/PyrusCommunis_BartlettDHv2.0.GO_GO_term.table.txt", head =T)

genus = "Pyrus"
species = "communisBartlettDH"  # Warning: the "_" is not allowed in the species name

makeOrgPackage(gene_info=gene_info,go=gene2go,
               version="0.1",
               maintainer = 'Y.N. <nieyuqi.cn@gmail.com>',
               author = 'Y.N. <nieyuqi.cn@gmail.com>',
               outputDir = output_dir,
               tax_id = "23211", # https://www.ncbi.nlm.nih.gov/taxonomy
               genus = genus,
               species = species,
               goTable="go")


## Genome 2
gene_info <- read.table("../output/s02.prepare_GO_KEGG/GWHBAOS00000000.GO_GID.table.txt", head =T)
gene2go <- read.table("../output/s02.prepare_GO_KEGG/GWHBAOS00000000.GO_GO_term.table.txt", head =T)

genus = "Pyrus"
species = "pyrifoliaCuiguan"

makeOrgPackage(gene_info=gene_info,go=gene2go,
               version="0.1",
               maintainer = 'Y.N. <nieyuqi.cn@gmail.com>',
               author = 'Y.N. <nieyuqi.cn@gmail.com>',
               outputDir = output_dir,
               tax_id = "3767", # https://www.ncbi.nlm.nih.gov/taxonomy
               genus = genus,
               species = species,
               goTable="go")
