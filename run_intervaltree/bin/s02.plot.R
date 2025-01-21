library(GenomicRanges)
library(ggplot2)

# Set working directory to the location of the script (if running interactively)
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}

rm(list=ls())

# Input -------------------------------------------------------------------

# Read input files
pre <-c ("GWHBAOS00000000", "PyrusCommunis_BartlettDHv2.0")

prefix <- pre[2]

prefix

gff <- read.table(paste0("../input/", prefix,".chr_list.gff"), header = FALSE, sep = "\t")
chr_lengths <- read.table(paste0("../input/", prefix,".genome.fasta.status.txt"), header = FALSE, sep = "\t")
accession_map <- read.table(paste0("../input/", prefix,".chrID_map.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)

output_file <- paste0("../output/s02.R.", prefix,".gene_density_pass.bed")

names(chr_lengths) <- c("chromosome", "length")
names(accession_map) <- c("accession", "chr_id")

# Sliding window parameters
window_size <- 1e6  # 1 Mb
step_size <- 1e5    # 100 kb

# Define cut-off values for different prefixes
cuto_offs <- data.frame(
  prefix = pre,
  value = c(30, 50)
)


# Set the cut-off value for the selected prefix
cut_off <- cuto_offs$value[cuto_offs$prefix == prefix]  # Extract the correct cut-off value based on the prefix

# Check the cut_off value
cut_off

# Calculating -------------------------------------------------------------

# Deal with the input files
genes <- gff[gff$V3 == "gene", ]

# Merge the ChrID map and set chr_id as a factor with correct levels
genes <- merge(genes, accession_map, by.x = "V1", by.y = "accession", all.x = TRUE)
genes$chr_id <- factor(genes$chr_id, levels = paste0("Chr", seq(1, 17)))

# Sort the genes dataframe by chr_id
genes <- genes[order(genes$chr_id), ]

# Set the chromosome factor levels correctly for chr_lengths as well
chr_lengths <- merge(chr_lengths, accession_map, by.x = "chromosome", by.y = "accession", all.x = TRUE)
chr_lengths$chr_id <- factor(chr_lengths$chr_id, levels = paste0("Chr", seq(1, 17)))

gr <- GRanges(seqnames = genes$chr_id, ranges = IRanges(start = genes$V4, end = genes$V5))
chr_length_vector <- setNames(chr_lengths$length, chr_lengths$chr_id)

# Sort chr_length_vector by level (Chr1, Chr2, ..., Chr17)
chr_length_vector <- chr_length_vector[levels(chr_lengths$chr_id)]

chr_length_vector

# Calculate offsets for global position
chr_offsets <- cumsum(c(0, head(chr_length_vector, -1)))
names(chr_offsets) <- names(chr_length_vector)



# Function to create sliding windows
create_sliding_windows <- function(chrom, chr_length, window_size, step_size) {
  starts <- seq(1, chr_length - window_size + 1, by = step_size)
  ends <- starts + window_size - 1
  GRanges(seqnames = chrom, ranges = IRanges(start = starts, end = ends))
}

# Create sliding windows
binned <- do.call(
  c, 
  lapply(names(chr_length_vector), function(chr) {
    create_sliding_windows(chr, chr_length_vector[chr], window_size, step_size)
  })
)

# Calculate coverage (overlaps with genes)
coverage <- countOverlaps(binned, gr)

# Create coverage data frame
coverage_df <- data.frame(
  start = start(binned),
  end = end(binned),
  coverage = as.numeric(coverage),
  chr = seqnames(binned)
)

# Add global start position
coverage_df$global_start <- coverage_df$start + chr_offsets[coverage_df$chr]


#  Plot -------------------------------------------------------------------


# Convert 'chr' to a factor with correct levels before plotting
coverage_df$chr <- factor(coverage_df$chr, levels = paste0("Chr", seq(1, 17)))


# Plot coverage with ChrID on the x-axis
ggplot(coverage_df) +
  geom_line(aes(x = global_start / 1e6, y = coverage, group = chr, color = chr)) +
  geom_hline(yintercept = cut_off, linetype = "dashed", color = "red") +  # Add cutoff line
  geom_rect(aes(xmin=0, xmax=sum(chr_lengths$length) / 1e6, ymin = 0, ymax = cut_off), 
           fill = "grey", alpha = 0.01) +  # Add grey rectangle from 0 to cut_off
  theme_minimal() +
  labs(
    title = paste0("Gene Coverage with Masked Regions (",prefix,")"),
    x = "Genome Position (Mb)",
    y = "Gene Count",
    color = "Chromosome"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )






# Write out a table -------------------------------------------------------

# Disable scientific notation
options(scipen = 999)

# Merge the ChrID map
coverage_df <- merge(coverage_df, accession_map, by.x = "chr", by.y = "chr_id", all.x = TRUE)

# Mask regions where coverage is less than the cutoff
coverage_df$masked <- coverage_df$coverage < cut_off  # mark low gene density window

# GRanges for masked regions
masked_gr <- GRanges(
  seqnames = coverage_df$chr[coverage_df$masked],
  ranges = IRanges(
    start = coverage_df$start[coverage_df$masked],
    end = coverage_df$end[coverage_df$masked]
  )
)

# Reduce the masked regions
reduced_masked_gr <- reduce(masked_gr)

# Calculate total masked length
masked_regions <- sum(width(reduced_masked_gr))  # the masked length after mask
genome_total_length <- sum(chr_length_vector)    # the total genome size
masked_percentage <- (masked_regions / genome_total_length) * 100

cat("The masked length:", masked_regions, "bp\n")
cat("The percentage of total genome size:", round(masked_percentage, 2), "%\n")

# Filter data for passing coverage
coverage_pass_df <- coverage_df[!coverage_df$masked, ]


# Create a new data frame for writing out with required columns
write_out <- data.frame(
  Chr = coverage_pass_df$accession,
  Start = coverage_pass_df$start - 1,  # Subtract 1 to get zero-based positions
  End = coverage_pass_df$end,
  Gene_count = coverage_pass_df$coverage
)

# Save the data frame with the pass/fail information to a CSV file
write.table(write_out, file = output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

