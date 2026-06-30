#!/usr/bin/env Rscript

## ================================================================
## Permutation test with regioneR — example using your test files
## ================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(123)

# 0) install + load
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(regioneR)
  library(GenomicRanges)
})

#TODO: add all genes

## 1) load region sets (your fictive test files)
special <- toGRanges("../data/PCOM/genes/positive_selection_summary_table.West.tsv.2kbp_upstream.bed")
altered <- toGRanges("../output/PCOM/MEGAnE_summary_PCOM.tsv.comm_Dessert.bed") # there are pop-specific TIPs
control <- toGRanges("../output/PCOM/MEGAnE_summary_PCOM.tsv.ALL.bed") # all TIPs


## genome from chrom.sizes
cs <- read.table(
  "../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat.sizes",
  sep = "\t", header = FALSE,
  col.names = c("chr", "size"), stringsAsFactors = FALSE
)
genome <- toGRanges(cs)
seqlengths(genome) <- end(genome)

## 2) permutation test: special vs altered
#Teste si les gènes sous sélection (special) chevauchent plus souvent les positions “altérées” (altered) que prévu au hasard.
pt1 <- permTest(
  A = special, B = altered,
  ntimes = 1000,
  alternative = "greater",
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions, #circularRandomizeRegions,
  genome = genome
)

## 3) permutation test: special vs control
#Même test mais avec ton échantillon témoin.
#Si le premier test est significatif et pas le second → signal réel.
pt2 <- permTest(
  A = special, B = control,
  ntimes = 1000,
  alternative = "greater",
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions,
  genome = genome
)

## 4) show summaries
cat("=== Genes under positive selection vs Dessert-specific TIPs ===\n")
print(pt1)
cat("=== Genes under positive selection vs Dessert-specific TIPs, reversed ===\n")
print(pt2)

## 5) plots
dir.create("../Figures/PCOM/regioneR", showWarnings = FALSE, recursive = TRUE)

# pt1 (Dessert)
png("../Figures/PCOM/regioneR/pear_permtests_Dessert_vs_ALL.png",
    width = 7, height = 6, units = "in", res = 300)
plot(pt1, main = "Permutation test: special vs Dessert")
dev.off()

# pt2 (Perry)
png("../Figures/PCOM/regioneR/pear_permtests_ALL_vs_Dessert.png",
    width = 7, height = 6, units = "in", res = 300)
plot(pt2, main = "Permutation test: special vs Perry")
dev.off()
