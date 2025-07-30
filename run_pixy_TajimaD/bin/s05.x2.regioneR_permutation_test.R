#!/usr/bin/env Rscript

###########################################################
# Permutation test with regioneR for genes under selection
###########################################################

#
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript this_script.R selected_genes.gff all_genes.bed genome_file.txt [output_prefix]")
}

selected_genes_file <- args[1]
all_genes_file <- args[2]
genome_file <- args[3]
output_dir <- args[4]

output_prefix <- tools::file_path_sans_ext(basename(selected_genes_file))

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rtracklayer)
  library(regioneR)
  library(GenomicRanges)
})

# Read input files
message("Input files:\n")
message(" - Selected genes:", selected_genes_file, "\n")
message(" - All genes:     ", all_genes_file, "\n")
message(" - Genome sizes:  ", genome_file, "\n")
message(" - Output to:     ", output_dir, "/", output_prefix, "\n\n")


if (endsWith(selected_genes_file, ".gff")) {
  message("Reading selected genes from GFF file...")
  selected_genes <- import(selected_genes_file)
} else if (endsWith(selected_genes_file, ".bed")) {
  message("Reading selected genes from BED file...")
  selected_genes <- toGRanges(selected_genes_file)
} else {
  stop("Unsupported file format for selected genes. Use GFF or BED.")
}

if (endsWith(all_genes_file, ".gff")) {
  message("Reading all genes from GFF file...")
  all_genes <- import(all_genes_file)
} else if (endsWith(all_genes_file, ".bed")) {
  message("Reading all genes from BED file...")
  all_genes <- toGRanges(all_genes_file)
} else {
  stop("Unsupported file format for all genes. Use GFF or BED.")
}

custom_genome  <- read.table(genome_file, header = FALSE, sep = "\t")
genome_def     <- toGRanges(custom_genome)



# Run permutation test
message("Running permutation test...")
pt <- permTest(
  A = selected_genes,
  B = all_genes,
  randomize.function = randomizeRegions, # Define randomization method
  genome = genome_def,
  evaluate.function = numOverlaps,
  ntimes = 1000,
  verbose = TRUE
)

# Pritn summary and p-value
message("Permutation test summary:")
print(pt)

pval <- pt$numOverlaps$pval
message("Permutation test p-value: ", pval)

# Save results
writeLines(paste0(pval), paste0(output_dir, "/", output_prefix, ".pvalue.txt"))

pdf(paste0(output_dir, "/", output_prefix, ".pdf"))
plot(pt)
dev.off()
