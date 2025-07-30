library(ggplot2)
library(tidyr)
library(dplyr)
library(rtracklayer)


# Set working directory to the location of the script and clean up
if (interactive()) {
	script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
	setwd(script_path)
	rm(list=ls())
}

###########################################################
# Permutation test with regioneR for genes under selection
###########################################################

# 1. Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("regioneR", quietly = TRUE))
    BiocManager::install("regioneR")

library(regioneR)
library(GenomicRanges)

###########################################################
# 2. Input data
###########################################################

# Replace these with your actual files:
# - genes_under_selection.bed: BED file with selected genes
# - all_genes.bed: BED file with all annotated genes

# use import() for gff files
# and toGRanges() for BED files
selected_genes <- import("../output/s04.catch_outliers/Pixy_TajimaD.0_0050outliers.betu.merged.bed.genes.gff")
#all_genes <- import("../input/GWHBAOS00000000.chr_list.gene.slop1kb.gff")
all_genes <- toGRanges("/shared/ifbstor1/projects/pear_snp3/pear/run_SIFT4G/input/PPY.cds_masked.bed")

# Genome information: you can create a genome size file like:
# chr1   50000000
# chr2   40000000
# ...
# Save as custom_genome.sizes
custom_genome <- read.table("../input/karyotype.txt", header = FALSE, sep = "\t")
genome_def <- toGRanges(custom_genome)

###########################################################
# 3. Define the randomization strategy
###########################################################

# Randomization will shuffle selected gene locations while preserving:
# - region lengths
# - avoiding genome gaps if defined
randomization_fun <- randomizeRegions

###########################################################
# 4. Perform the permutation test
###########################################################

pt <- permTest(
  A = selected_genes,            # Regions of interest (genes under selection)
  B = all_genes,                 # Background gene set
  randomize.function = randomization_fun,
  genome = genome_def,           # Genome definition (chrom sizes)
  evaluate.function = numOverlaps, # Statistic: count overlaps
  ntimes = 1000,                 # Number of permutations
  verbose = TRUE
)

###########################################################
# 5. Inspect results
###########################################################

# Print summary
print(pt)

# Plot the permutation distribution
plot(pt)

# Extract p-value
cat("Permutation test p-value:", pt$numOverlaps$pval, "\n")

###########################################################
# 6. Interpretation
###########################################################
# - If p-value < 0.05: the overlap between selected genes and background
#   is significantly different from random expectation.
# - If p-value > 0.05: cannot reject the null hypothesis that
#   the observed overlap is due to chance.