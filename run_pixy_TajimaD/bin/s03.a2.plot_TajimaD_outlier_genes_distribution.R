# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)


# Set working directory to the location of the script and clean up
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}

rm(list = ls())

# Read real data ----------------------------------------------------------
cutoff_str <- "0_0040"
cutoff_val <- 0.004
cutoff_pct <- '0.4%'

# Dir for real Tajima D output files
base_dir1 <- "../../run_pixy_TajimaD"
base_dir2 <- "../../run_pixy_TajimaD_refer_on_communis"
base_dirs <- c(base_dir1, base_dir2)
folder <- c("Oriental", "Occidental")

all_counts_long <- list()
all_pval_df <- list()

for (inx in c(1, 2)) {

masked_dir  <- file.path(base_dirs[inx], "output", "s04.catch_outliers")
masked_pattern  <- "Pixy_TajimaD\\.0_\\d{4}outliers\\.(.*?)\\.merged\\.bed\\.genes\\.txt$"

masked_files  <- list.files(masked_dir,  pattern = masked_pattern, full.names = TRUE)

# Function to extract counts from files
extract_counts <- function(files, type, pattern) {
  tibble(File = files) %>%
    mutate(
      Filename   = basename(File),
      Cutoff_raw = str_extract(Filename, "0_\\d{4}"),
      Cutoff     = as.numeric(sub("0_", "0.", Cutoff_raw)),
      Population = str_match(Filename, pattern)[, 2],
      Count      = map_int(File, ~ length(readLines(.x))),
      Type       = type
    ) %>%
    select(Population, Cutoff_raw, Cutoff, Type, Count)
}

# Extract counts from raw and masked files
masked_df  <- extract_counts(masked_files,  "Masked",    masked_pattern)

# Combine and arrange counts data
counts_long <- bind_rows(masked_df) %>%
  arrange(Population, Cutoff, Type) %>%
  mutate(
    linewidth = ifelse(Type == "Masked", 1, 0.4),
    Group = folder[inx]
  )

# Define p value dirs
pval_dir <- file.path(base_dirs[inx], "output", "s05.regioneR_permutation_test")
pval_files <- list.files(pval_dir, 
    pattern = "^Pixy_TajimaD\\.0_\\d{4}outliers\\..*?\\.pvalue\\.txt$", 
    full.names = TRUE)

# get p value data
pval_df <- tibble(File = pval_files) %>%
  mutate(
    Filename = basename(File),
    Cutoff_raw = str_extract(Filename, "0_\\d{4}"),
    Cutoff = as.numeric(sub("0_", "0.", Cutoff_raw)),
    Population = str_replace(Filename, 
      "^Pixy_TajimaD\\.0_\\d{4}outliers\\.(.*?)\\.merged\\.bed\\.genes\\.pvalue\\.txt$", "\\1"),
    Pvalue = map_dbl(File, ~ as.numeric(readLines(.x, n = 1))),
    P_label = paste0("p=", formatC(Pvalue, format = "f", digits = 4)),
    Group = folder[inx]
  )

# store results from West and East in lists
all_counts_long[[inx]] <- counts_long
all_pval_df[[inx]] <- pval_df

# Create labels for the plot
label_df <- counts_long %>%
  filter(Cutoff_raw == cutoff_str)

label_df <- label_df %>%
  left_join(pval_df, by = c("Population", "Cutoff_raw")) %>%
  mutate(
    x_pos = Cutoff.x + 0.0003,
    y_pos = log10(Count) - 0.1,
    y_ppos = Pvalue + 0.004,
    label_text = paste0(Count, " genes  ", P_label),
    gene_label = paste0(Count, " genes at cutoff ", cutoff_pct),
    p_label = paste0("p=", formatC(Pvalue, format = "f", digits = 4),
                     " at cutoff ", cutoff_pct)
  )

head(label_df)

# Load color configuration for pouplation ------------------------------
# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
ori_order_color <- setNames(color_uniq$Color, color_uniq$Population)

pops_in_data <- unique(counts_long$Population)
my_colors <- ori_order_color[pops_in_data]


# ----------------------------------
#  plot outlier gene counts ~ cutoff
# ----------------------------------

# Combined plot
p1 <- ggplot(counts_long, aes(x = Cutoff, y = log10(Count),
                              group = interaction(Population, Type),
                              color = Population)) +
  geom_line(aes(linewidth = linewidth), linetype = "solid", alpha = 0.6) +
  geom_vline(xintercept = cutoff_val, linetype = "dashed", color= "black",
            size = 0.4, alpha = 0.8)+
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values = my_colors) +
  scale_linewidth_identity() +
  labs(title = paste0("Outlier Gene Counts at Different Tajima's D Cutoffs (", folder[inx], ")"),
       x = "Cutoff", y = "log10(Number of Genes)") +
  theme_minimal(base_size = 6) +
  theme(
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Faceted plot per population
p2 <- ggplot(counts_long, aes(x = Cutoff, y = log10(Count),
                              group = interaction(Population, Type),
                              color = Population)) +
  geom_line(aes(linewidth = linewidth), linetype = "solid", alpha = 0.6) +
  geom_point(size = 1, shape = 16) +

  geom_vline(xintercept = cutoff_val, linetype = "dashed", color= "black",
             size = 0.4, alpha = 0.8)+
  geom_text(data = label_df, aes(x = x_pos, y = y_pos, label = gene_label),
            size = 1.6, vjust = 1, hjust=0, show.legend=FALSE) +

  facet_wrap(~Population, nrow = 4, ncol = 4) +
  scale_color_manual(values = my_colors) +
  scale_linewidth_identity() +
  labs(title = "Outlier Gene Counts per Population",
       x = "Cutoff", y = "log10(Number of Genes)") +
  theme_minimal(base_size = 6) +
  theme(
    # panel.grid = element_blank(),
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine and save
final_plot <- p1 / p2 + plot_layout(heights = c(1, 2))

final_plot


ggsave(
  filename = paste0("s03.a2.TajimaD_GeneCount_Distribution_", folder[inx], ".png"),
  plot = final_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600
)

ggsave(
  filename = paste0("s03.a2.TajimaD_GeneCount_Distribution_", folder[inx], ".tif"),
  plot = final_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

# -----------------------
#  plot pvalue ~ cutoff 
# -----------------------

# global p value
pval_p1 <- ggplot(pval_df, aes(x = Cutoff, y = Pvalue, color = Population, group = Population)) +
  geom_line(linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 1, shape = 16) +
  geom_vline(xintercept = cutoff_val, linetype = "dashed", color = "black",
             size = 0.4, alpha = 0.8) +
  scale_color_manual(values = my_colors) +
  labs(title = paste0("Permutation Test P-values at Different Tajima's D Cutoffs (", folder[inx], ")"),
       x = "Cutoff", y = "P-value") +
  theme_minimal(base_size = 6) +
  theme(
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# per population facet
pval_p2 <- ggplot(pval_df, aes(x = Cutoff, y = Pvalue, color = Population, group = Population)) +
  geom_line(linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 1, shape = 16) +
  geom_vline(xintercept = cutoff_val, linetype = "dashed", color = "black",
             size = 0.4, alpha = 0.8) +
  geom_text(data = label_df, aes(x = x_pos, y = y_ppos, label = p_label),
            size = 1.6, vjust = 1, hjust=0, show.legend=FALSE) +
   geom_text(data = label_df, aes(x = x_pos, y = 0.03, label = ""), show.legend=FALSE) +

  facet_wrap(~Population, nrow = 4, ncol = 4) +
  scale_color_manual(values = my_colors) +
  labs(title = "Permutation Test P-values per Population",
       x = "Cutoff", y = "P-value") +
  theme_minimal(base_size = 6) +
  theme(
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# final pvalue 
final_pval_plot <- pval_p1 / pval_p2 + plot_layout(heights = c(1, 2))
final_pval_plot

# Save pvalue 
ggsave(
  filename = paste0("s03.a2.TajimaD_Pvalue_Distribution_", folder[inx], ".png"),
  plot = final_pval_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600
)

ggsave(
  filename = paste0("s03.a2.TajimaD_Pvalue_Distribution_", folder[inx], ".tif"),
  plot = final_pval_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)
}


# -------------------
# Save tables
# -------------------
all_counts_long_df <- bind_rows(all_counts_long, .id = "Group")
all_pval_df_df <- bind_rows(all_pval_df, .id = "Group")

counts_wide <- all_counts_long_df %>%
  pivot_wider(
    names_from = Cutoff,
    values_from = Count,
    names_sort = TRUE,
    names_prefix = "TajimaD_cutoff_"
  )
write.table(counts_wide, "s03.a2.TajimaD_GeneCount_Table_All.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

pval_wide <- all_pval_df_df %>%
  pivot_wider(
    names_from = Cutoff,
    values_from = Pvalue, 
    names_sort = TRUE,
    names_prefix = "TajimaD_cutoff_"
  )
write.table(pval_wide, "s03.a2.TajimaD_Pvalue_WideTable_All.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
