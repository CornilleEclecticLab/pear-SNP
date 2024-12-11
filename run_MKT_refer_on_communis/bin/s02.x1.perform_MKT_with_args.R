#!/usr/bin/env Rscript

# Reference:
# https://eacooper400.github.io/gen8900/exercises/mk.html

# Load required libraries
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(Rsamtools)
  library(tidyverse)
  library(data.table)
})

# Function to perform McDonald-Kreitman Test
perform_mkt_analysis <- function(
    vcf_file,                 ## Path to VCF file
    gff_file,                 ## Path to GFF file
    fasta_file,               ## Path to FASTA file
    sample_file,              ## Path to sample file
    outgroup_pop = NULL,      ## Outgroup population 
    output_dir = "../output"  ## Directory for output files
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load sample information
  pop_samples <- fread(sample_file, header = FALSE, col.names=c("sample", "population"))
  populations <- unique(pop_samples$population)
  
  # Set outgroup population (if not specified, use first population)
  if (is.null(outgroup_pop)) {
    outgroup_pop <- populations[1]
  }
  
  # Open reference files
  fa <- open(FaFile(fasta_file))
  txdb <- makeTxDbFromGFF(gff_file)
  
  # Read VCF file
  vcf.object <- readVcf(file = vcf_file)
  
  # Variant annotation
  effects <- predictCoding(vcf.object, txdb, fa)
  
  # Match VCF entries with annotated effects
  id.from.effects <- names(ranges(effects))
  id.from.vcf <- rownames(info(vcf.object))
  m <- match(id.from.vcf, id.from.effects)
  
  # Assign gene IDs based on matching
  gene_ids <- rep(".", nrow(vcf.object))
  gene_ids[!is.na(m)] <- effects$GENEID[m[!is.na(m)]]
  rownames(vcf.object) <- gene_ids
  
  # Annotate VCF with the CONSEQUENCE field
  info(vcf.object)$CSQ <- effects$CONSEQUENCE[m]
  
  # Filter VCF to keep only entries with matching annotations
  vcf.annotated <- vcf.object[!is.na(m)]
  rm(vcf.object)
  
  # Read VCF data into a dataframe
  read.vcf <- function(file, special.char="##", ...) {
    my.search.term = paste0(special.char, ".*")
    all.lines = readLines(file)
    clean.lines = gsub(my.search.term, "", all.lines)
    clean.lines = gsub("#CHROM", "CHROM", clean.lines)
    read.table(..., text=paste(clean.lines, collapse="\n"))
  }
  
  vcf_filebase <- tools::file_path_sans_ext(vcf_file)
  middle_annotated_file <- file.path(output_dir, paste0(vcf_filebase, "_MKT_Ann.vcf"))
  my.data <- read.vcf(writeVcf(vcf.annotated, middle_annotated_file),
                      header = TRUE, stringsAsFactors = FALSE)
  my.data$INFO <- gsub(".*;CSQ=", "", my.data$INFO)
  rm(vcf.annotated)
  
  # Some functions to deal with VCF 
  ## Credit: https://gist.github.com/eacooper400/c954c757db02eb6c4ccfae1aa090658c
  allele.freq <- function(genotypeCounts) {
    n = sum(genotypeCounts) - genotypeCounts["NN"]
    p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
    q = 1 - p
    freqs = c(p, q)
    names(freqs) = c("p", "q")
    return(freqs)
  }
  
  count.genotypes <- function(genotypes) {
    genotypes = gsub("(\\||/)", "", genotypes) 
    gen.patterns = c("00", "01", "10", "11", "..") 
    my.counts = table(factor(genotypes, levels=gen.patterns)) 
    final.counts = c(my.counts[1], (my.counts[2] + my.counts[3]), my.counts[4:5]) 
    names(final.counts) = c("AA", "Aa", "aa", "NN") 
    return(final.counts)
  }
  
  get.field <- function(samples, format, fieldName) {
    x = strsplit(samples, split=":")
    fields = unlist(strsplit(format, split=":")) 
    i = which(fields == fieldName)
    if (!(fieldName %in% fields)) stop('fieldName not found in format fields') 
    return(sapply(x, `[[`, i)) 
  }
  
  # Prepare data for analysis
  new <- my.data[, c(3, 8)]
  
  # Identify the genotypes for outgroup
  new[[paste0(outgroup_pop, "_refAF")]] <- NA
  pop_samples_list <- pop_samples$sample[pop_samples$population == outgroup_pop]
  pop.data <- my.data[, c(1:9, which(colnames(my.data) %in% pop_samples_list))]
  
  for (i in 1:nrow(pop.data)) {
    my.row <- as.vector(pop.data[i, ], mode="character")
    genotypes <- get.field(my.row[10:length(my.row)], my.row[9], "GT")
    all.counts <- count.genotypes(genotypes)
    ref.freq <- allele.freq(all.counts)["p"]
    new[[paste0(outgroup_pop, "_refAF")]][i] <- ref.freq
  }
  # No need classification for outgroup population
  # new[[paste0(outgroup_pop, "_SiteClass")]] <- ifelse(
  #  abs(new[[paste0(outgroup_pop, "_refAF")]]) == 1 | abs(new[[paste0(outgroup_pop, "_refAF")]]) == 0,
  #  "Fixed", 
  #  "Polymorphic"
  # )
  
  # Identify genotypes for each population
  for (pop in populations) {
    if (pop == outgroup_pop) next
    new[[paste0(pop, "_refAF")]] <- NA
    pop_samples_list <- pop_samples$sample[pop_samples$population == pop]
    pop.data <- my.data[, c(1:9, which(colnames(my.data) %in% pop_samples_list))]
    
    for (i in 1:nrow(pop.data)) {
      my.row <- as.vector(pop.data[i, ], mode="character")
      genotypes <- get.field(my.row[10:length(my.row)], my.row[9], "GT")
      all.counts <- count.genotypes(genotypes)
      ref.freq <- allele.freq(all.counts)["p"]
      new[[paste0(pop, "_refAF")]][i] <- ref.freq
    }
    
    new[[paste0(pop, "_SiteClass")]] <- ifelse(
      abs(new[[paste0(pop, "_refAF")]] - new[[paste0(outgroup_pop, "_refAF")]]) == 1, 
      "Fixed", 
      "Polymorphic"
    )
  }
  
  # Calculate Results: N:S for Between and Within each populations and outgroup pop
  geneList <- unique(my.data$ID)
  my.results.allpop <- list()
  
  for (pop in populations) {
    if (pop == outgroup_pop) next  # Skip outgroup population
    
    my.results <- data.frame(
      Genes = geneList, 
      Fn = 0,  ## Fixed nonsynonymous differences
      Fs = 0,  ## Fixed synonymous differences
      Pn = 0,  ## Polymorphic nonsynonymous sites
      Ps = 0,  ## Polymorphic synonymous sites
      stringsAsFactors = FALSE
    )
    
    for (i in 1:length(geneList)) {
      gene <- new[new$ID == geneList[i], ]
      if(nrow(gene) > 0) {
        site_class_col <- paste0(pop, "_SiteClass")
        
        fixed_syn <- gene[gene$INFO == "synonymous" & 
                            gene[[site_class_col]] == "Fixed", ]
        fixed_non <- gene[gene$INFO != "synonymous" &  ## Do not use =="nonsynonymous"
                                                       ## or miss other types, e.g. frameshift
                            gene[[site_class_col]] == "Fixed", ]
        
        poly_syn <- gene[gene$INFO == "synonymous" & 
                           gene[[site_class_col]] == "Polymorphic", ]
        poly_non <- gene[gene$INFO != "synonymous" & 
                           gene[[site_class_col]] == "Polymorphic", ]
        
        my.results$Fn[i] <- nrow(fixed_non)
        my.results$Fs[i] <- nrow(fixed_syn)
        my.results$Pn[i] <- nrow(poly_non)
        my.results$Ps[i] <- nrow(poly_syn)
      }
    }
    
    # Perform Fisher's exact test and calculate additional metrics
    my.results$pValues <- apply(my.results[, 2:5], 1, function(x) {
      contingency_matrix <- matrix(x, nrow = 2, byrow = TRUE)
      tryCatch(fisher.test(contingency_matrix)$p.value, error = function(e) NA)
    })
    
    # Calculate neutrality indices
    my.results$BetweenN_S <- ifelse(my.results$Fs > 0, 
                                    my.results$Fn / my.results$Fs, 
                                    NA)
    my.results$WithinN_S <- ifelse(my.results$Ps > 0, 
                                   my.results$Pn / my.results$Ps, 
                                   NA)
    my.results$NI <- ifelse(my.results$Fs > 0 & my.results$Ps > 0, 
                            (my.results$Pn / my.results$Ps) / (my.results$Fn / my.results$Fs), 
                            NA)
    
    # Store results for each population
    my.results.allpop[[pop]] <- my.results
  }
  
  # Write results in wide format
  write_mkt_results_wide <- function(results_list, output_dir, vcf_file) {
    # Get unique genes across all populations
    all_genes <- unique(unlist(lapply(results_list, function(x) x$Genes)))
    
    # Create a wide-format data frame
    wide_results <- data.frame(Gene = all_genes)
    
    # Iterate through populations to add columns
    for (pop in names(results_list)) {
      # Get results for this population
      pop_results <- results_list[[pop]]
      
      # Create columns for this population's p-values and NI
      p_value_col <- paste0(pop, "_Pvalues")
      ni_col <- paste0(pop, "_NI")
      
      # Merge results with wide_results
      wide_results <- merge(wide_results, 
                            pop_results[, c("Genes", "pValues", "NI")], 
                            by.x = "Gene", 
                            by.y = "Genes", 
                            all.x = TRUE)
      
      # Rename columns
      names(wide_results)[names(wide_results) == "pValues"] <- p_value_col
      names(wide_results)[names(wide_results) == "NI"] <- ni_col
    }
    
    # Create filename for wide-format results
    vcf_filename <- basename(vcf_file)
    vcf_basename <- tools::file_path_sans_ext(vcf_filename)
    filename <- file.path(output_dir, paste0(vcf_basename, "_MKT_wide_results.csv"))
    
    # Write the wide-format results to CSV
    write.csv(wide_results, filename, row.names = FALSE, quote = FALSE)
    
    cat("McDonald-Kreitman Test results in wide format written to", filename, "\n")
  }
  
  # Call function to write results
  write_mkt_results_wide(my.results.allpop, output_dir, vcf_file)
  
  # Return results for potential further analysis
  return(my.results.allpop)
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) < 5) {
  cat("Usage: Rscript MKT_analysis.R <vcf_file> <gff_file> <fasta_file> <sample_file> <outgroup_pop>\n")
  quit(status = 1)
}

# Assign arguments to variables
vcf_file <- args[1]
gff_file <- args[2]
fasta_file <- args[3]
sample_file <- args[4]
outgroup_pop <- args[5]

# Validate input files
input_files <- c(vcf_file, gff_file, fasta_file, sample_file)
for (file in input_files) {
  if (!file.exists(file)) {
    cat("Error: File does not exist -", file, "\n")
    quit(status = 1)
  }
}

# Create log file
log_file <- file.path("../output", paste0("MKT_analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
dir.create("../output", showWarnings = FALSE)
sink(log_file, append = TRUE)

# Start timing
start_time <- Sys.time()
cat("McDonald-Kreitman Test Analysis Started\n")
cat("Date:", as.character(start_time), "\n")
cat("Input Files:\n")
cat("  VCF File:", vcf_file, "\n")
cat("  GFF File:", gff_file, "\n")
cat("  FASTA File:", fasta_file, "\n")
cat("  Sample File:", sample_file, "\n")
cat("  Outgroup Population:", outgroup_pop, "\n")

# Run analysis
tryCatch({
  results <- perform_mkt_analysis(
    vcf_file = vcf_file,
    gff_file = gff_file,
    fasta_file = fasta_file,
    sample_file = sample_file,
    outgroup_pop = outgroup_pop,
    output_dir = "../output"
  )
  
  # End timing
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  cat("\nAnalysis Completed Successfully\n")
  cat("Total Runtime:", as.numeric(duration), "minutes\n")
  
  # Close log file
  sink()
  
  # Exit with success status
  quit(status = 0)
}, error = function(e) {
  # Error handling
  cat("Error in McDonald-Kreitman Test Analysis:\n")
  cat(conditionMessage(e), "\n")
  
  # Close log file
  sink()
  
  # Exit with error status
  quit(status = 1)
})