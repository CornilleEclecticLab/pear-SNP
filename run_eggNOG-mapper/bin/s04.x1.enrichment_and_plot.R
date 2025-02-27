# By Xilong CHEN
# Created on: 2023-01-25
# Contact: chen_xilong@outlook.com
# 
# Modified by Yuqi NIE on 2025-02-11

# Set Working Directory to source file folder
rm(list = ls())

# Load required libraries
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(AnnotationForge)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

database <- args[1]  # "comm" or "pyri"
gene_list_f <- args[2]
output_dir <- args[3]


# Validate database selection
if (database == "comm") {
    db <- "org.PcommunisBartlettDH.eg.db"
    kegg_db_prefix <- "PyrusCommunis_BartlettDHv2.0"
} else if (database == "pyri") {
    db <- "org.PpyrifoliaCuiguan.eg.db"
    kegg_db_prefix <- "GWHBAOS00000000"
} else {
    stop("Error: Please input a valid database name: 'comm' or 'pyri'.")
}

# Load GO and KEGG database
go_base_path <- file.path("../output/s03.build_GO_KEGG_database", db)
install.packages(go_base_path, repos = NULL, type = "source")
library(as.character(db), character.only = TRUE)

# Read gene list
gene_df <- read.table(gene_list_f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes <- gene_df$V1

# GO enrichment analysis
GO <- enrichGO(
    genes,
    OrgDb = db,
    keyType = "GID",
    ont = "ALL",
    pvalueCutoff = 0.05,
    minGSSize = 1,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = FALSE
)

# GO enrichment plot
p_GO <- dotplot(GO, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scales = "free")

# Read KEGG pathway mapping files
pathway2gene <- read.table(
    file.path("../output/s02.prepare_GO_KEGG", paste0(kegg_db_prefix, ".KEGG_TERM2GENE.table.txt")),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

pathway2name <- read.table(
    "../output/s02.prepare_GO_KEGG/KEGG_TERM2NAME.table.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# KEGG enrichment analysis
KEGG <- enricher(
    as.vector(genes),
    TERM2GENE = pathway2gene,
    TERM2NAME = pathway2name,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 1
)

# KEGG enrichment plot
p_KEGG <- dotplot(KEGG)

# Save plots
basename_f <- basename(gene_list_f)
ggsave(file.path(output_dir, paste0(basename_f, ".GO_enrichment.pdf")), p_GO, width = 6, height = 10)
ggsave(file.path(output_dir, paste0(basename_f, ".KEGG_enrichment.pdf")), p_KEGG, width = 8, height =6)

# Save enrichment analysis results
write.table(as.data.frame(GO), file.path(output_dir, paste0(basename_f, ".GO_enrichment.txt")), 
            sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(KEGG), file.path(output_dir, paste0(basename_f, ".KEGG_enrichment.txt")), 
            sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
