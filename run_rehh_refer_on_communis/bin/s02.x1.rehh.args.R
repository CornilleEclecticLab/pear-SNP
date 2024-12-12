suppressPackageStartupMessages({
    library(rehh)
    library(data.table)
    library(dplyr)
    library(parallel)
    library(optparse)
})

# Parse command-line arguments
option_list <- list(
    make_option(c("-s", "--samples"), type = "character",
                help = "Path to the population samples file", metavar = "character"),
    make_option(c("-o", "--outgroup"), type = "character", 
                help = "Outgroup population name", metavar = "character"),
    make_option(c("-v", "--vcf_folder"), type = "character",
                help = "Folder containing VCF files", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Assign parsed arguments to variables
pop_samples_path <- opt$samples
outgroup_pop <- opt$outgroup
vcf_folder <- opt$vcf_folder


# Set working directory to the location of the script
# script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(script_path)
# rm(list = ls())

# Output directory
output_dir <- "../output/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Function to write results to TSV files
write_table <- function(data, file_name) {
    write.table(data, file = paste0(output_dir, file_name), 
                sep = "\t", quote = FALSE, row.names = TRUE)
}

# Function to save plots to PDF
save_plot <- function(file_name, plot_code) {
    pdf(file_name)
    plot_code()
    dev.off()
}

# Read the population files
pop_samples <- read.table(pop_samples_path, 
                          header = FALSE, sep = '\t', 
                          col.names = c("sample", "population"))
populations <- unique(pop_samples$population)
rm(pop_samples)


# Ensure the outgroup is included in total populations
if (!(outgroup_pop %in% populations)) {
    populations <- c(populations, outgroup_pop)
}

# List VCF files
vcf_files_paths <- list.files(vcf_folder, pattern = "masked.vcf$", full.names = TRUE)

# Extract chromosome number and population from file names
extract_chromosome_number <- function(file_path) {
    sub("(.*)\\.(.*)\\.(.*)\\.masked\\.vcf", "\\2", basename(file_path))
}

extract_population_name <- function(file_path) {
    sub("(.*)\\.(.*)\\.(.*)\\.masked\\.vcf", "\\3", basename(file_path))
}

# Data frame of VCF files paths, chromosome ID, and populations
vcf_files <- vcf_files_paths %>%
    lapply(function(file) {
        data.frame(
            File = file,
            Chromosome = extract_chromosome_number(file),
            Population = extract_population_name(file)
        )
    }) %>%
    bind_rows() %>%
    filter(Population %in% populations)

# Ensure consistent VCF file counts across populations
vcf_n <- 0
for (pop in populations) {
    pop_vcf_count <- sum(vcf_files$Population == pop)
    if (pop_vcf_count == 0) {
        stop(paste("No VCF files found for population", pop))
    }
    if (vcf_n == 0) {
        vcf_n <- pop_vcf_count
    } else if (pop_vcf_count != vcf_n) {
        stop("Different number of VCF files for populations")
    }
}

chrs <- unique(vcf_files$Chromosome)



# Run rehh scans
wgscan_pops <- list()
for (pop in populations) {
    wgscan_pop <- NULL
    for (chr in chrs) {
        vcf <- pull(filter(vcf_files, 
                           Population == pop, Chromosome == chr), File)

        # Create internal haplotype representation
        hh <- data2haplohh(hap_file = vcf, polarize_vcf = FALSE, vcf_reader = "data.table")

        # Perform scan on a single chromosome
        chr_scan <- scan_hh(hh)

        # Combine chromosome-wise data frames
        wgscan_pop <- rbind(wgscan_pop, chr_scan)
    }
    write_table(wgscan_pop, paste0("s02.", pop, ".wgscan.tsv"))
    wgscan_pops[[pop]] <- wgscan_pop
}


# Calculate genome-wide iHS values
ihs_pops <- list()
for (pop in populations) {
    wgscan_pop <- wgscan_pops[[pop]]
    ihs_pop <- ihh2ihs(wgscan_pop, 
                       freqbin = 1) ## For unpolarized data, freqbin = 1
                                    ## Because the standardization of iHS
                                    ## should not be done in a frequency-dependent manner
    write_table(ihs_pop$ihs, paste0("s02.", pop, ".ihs.tsv"))
    ihs_pops[[pop]] <- ihs_pop
}

# Calculate Rsb values with the outgroup
rsb_pops <- list()
for (pop in setdiff(populations, outgroup_pop)) {
    rsb <- ines2rsb(scan_pop1 = wgscan_pops[[outgroup_pop]],
                    scan_pop2 = wgscan_pops[[pop]],
                    popname1 = outgroup_pop,
                    popname2 = pop,
                    include_freq = TRUE,
                    p.side = 'left')  ## Default is two-sided 
                                      ## Set "left" to identify strong extended homozygosity 
                                      ## in population pop2 relative to population pop1 (outgroup)
    write_table(rsb, paste0("s02.", pop, ".rsb.tsv"))
    rsb_pops[[pop]] <- rsb
}


# Calculate XP-EHH values with the outgroup
xpehh_pops <- list()
for (pop in setdiff(populations, outgroup_pop)) {
    xpehh <- ies2xpehh(scan_pop1 = wgscan_pops[[outgroup_pop]],
                       scan_pop2 = wgscan_pops[[pop]],
                       popname1 = outgroup_pop,
                       popname2 = pop)
    write_table(xpehh, paste0("s02.", pop, ".xpehh.tsv"))
    xpehh_pops[[pop]] <- xpehh
}


# Visualize Rsb vs. XP-EHH
for (pop in populations) {
    if (pop == outgroup_pop) next

    rsb <- rsb_pops[[pop]]
    xpehh <- xpehh_pops[[pop]]

    save_plot(paste0(output_dir, "s02.", pop, ".rsb_vs_xpehh_plot.pdf"), 
        function() {
            plot(rsb[, paste0("RSB_", outgroup_pop, "_", pop)],
                 xpehh[, paste0("XPEHH_", outgroup_pop, "_", pop)],
                 xlab = "Rsb",
                 ylab = "XP-EHH",
                 pch = ".",
                 xlim = c(-7.5, 7.5),
                 ylim = c(-7.5, 7.5))

            # Add dashed diagonal line
            abline(a = 0, b = 1, lty = 2)
    })
}

# Visualize iHS Distribution and Manhattan plots
for (pop in populations) {
    ihs_pop <- ihs_pops[[pop]]
    ihs <- ihs_pop$ihs$IHS

    save_plot(paste0(output_dir, "s02.", pop, ".ihs_distplot.pdf"), function() {
        distribplot(ihs, 
                 xlab = "iHS", 
                 main = paste0("Distribution of iHS values for ", pop))
    })

    save_plot(paste0(output_dir, "s02.", pop, ".ihs_qqplot.pdf"), function() {
        distribplot(ihs, 
                    xlab = "iHS", 
                    qqplot = TRUE,
                    main = paste0("Q-Q Plot of iHS values for ", pop))
    })

    save_plot(paste0(output_dir, "s02.", pop, ".ihs_manhattanplot.pdf"), function() {
        manhattanplot(ihs_pop, 
                      main = paste0("iHS values for ", pop))
    })

    save_plot(paste0(output_dir, "s02.", pop, ".ihs_manhattanplot_pval.pdf"), function() {
        manhattanplot(ihs_pop, 
                      pval = TRUE, 
                      threshold = 4, 
                      main = paste0("p-value of iHS (", pop, ")"))
    })
}
