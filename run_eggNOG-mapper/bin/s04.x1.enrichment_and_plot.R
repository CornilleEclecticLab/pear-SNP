# By Xilong CHEN
# Created on: 2023-01-25
# Contact: chen_xilong@outlook.com
# 
# Modified by Yuqi NIE on 2025-02-11

# Set Working Directory to source file folder
rm(list = ls())

# Load required libraries
suppressMessages({
    library(ggplot2)
    library(stringr)
    library(enrichplot)
    library(clusterProfiler)
    library(AnnotationForge)
})


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

database <- args[1]  # "comm" or "pyri"
gene_list_f <- args[2]
output_dir <- args[3]

message("Database selected: ", database)
message("Gene list file: ", gene_list_f)
message("Output directory: ", output_dir)


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

if (!requireNamespace(db, quietly = TRUE)) {
    install.packages(go_base_path, repos = NULL, type = "source")
}

library(as.character(db), character.only = TRUE)

# Read gene list
gene_df <- read.table(gene_list_f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes <- gene_df$V1

# Initialize variables to avoid "object not found" errors
GO <- NULL
go_simplified <- NULL
p_GO <- NULL
p_go_sim <- NULL
KEGG <- NULL
p_KEGG <- NULL

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
    pool = FALSE,
    readable = FALSE    # If readable is set to TRUE, the input gene IDs will be converted to gene symbols.
)

# Reduce redundancy in GO terms
if (is.null(GO) || nrow(as.data.frame(GO)) == 0) {
    message("No GO terms enriched.")
} else {
    message("GO enrichment analysis completed.")

    # Check ONTOLOGY column
    if (!"ONTOLOGY" %in% colnames(as.data.frame(GO))) {
        GO@result$ONTOLOGY <- "Unknown"
    }

    go_simplified <- simplify(GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
    if (!"ONTOLOGY" %in% colnames(as.data.frame(go_simplified))) {
        go_simplified@result$ONTOLOGY <- "Unknown"
    }

    # GO enrichment plots
    p_GO <- dotplot(GO, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scales = "free")

    p_go_sim <- dotplot(go_simplified, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scales = "free")
}


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
if (!is.null(KEGG) && nrow(as.data.frame(KEGG)) > 0) {
    p_KEGG <- dotplot(KEGG)
} else {
    message("No KEGG terms enriched.")
    p_KEGG <- NULL
}

# Save plots
basename_f <- basename(gene_list_f)
if (!is.null(p_GO)) {
    ggsave(file.path(output_dir, paste0(basename_f, ".GO_enrichment.pdf")), p_GO, width = 6, height = 10)
}
if (!is.null(p_go_sim)) {
    ggsave(file.path(output_dir, paste0(basename_f, ".GO_simplified_enrichment.pdf")), p_go_sim, width = 6, height = 10)
}
if (!is.null(p_KEGG)) {
    ggsave(file.path(output_dir, paste0(basename_f, ".KEGG_enrichment.pdf")), p_KEGG, width = 8, height = 6)
}


# Get population name from gene list file
parts <- strsplit(gene_list_f, "\\.")[[1]]
Population <- parts[length(parts) - 1]

# Add population name to results and save into tables
if (!is.null(GO) && nrow(as.data.frame(GO)) > 0) {
    GO_out <- cbind(Population = Population, as.data.frame(GO))
    write.table(GO_out, file.path(output_dir, paste0(basename_f, ".GO_enrichment.txt")), 
                sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
}

if (!is.null(go_simplified) && nrow(as.data.frame(go_simplified)) > 0) {
    GO_simple_out <- cbind(Population = Population, as.data.frame(go_simplified))
    write.table(GO_simple_out, file.path(output_dir, paste0(basename_f, ".GO_simplified_enrichment.txt")), 
                sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
}

if (!is.null(KEGG) && nrow(as.data.frame(KEGG)) > 0) {
    KEGG_out <- cbind(Population = Population, as.data.frame(KEGG))
    write.table(KEGG_out, file.path(output_dir, paste0(basename_f, ".KEGG_enrichment.txt")), 
                sep = "\t", na = "nan", quote = FALSE, row.names = FALSE)
}
