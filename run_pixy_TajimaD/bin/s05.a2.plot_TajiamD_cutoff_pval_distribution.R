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

# Dir for real Tajima D output files
base_dir1 <- "../../run_pixy_TajimaD"
base_dir2 <- "../../run_pixy_TajimaD_refer_on_communis"
base_dirs <- c(base_dir1, base_dir2)
folder <- c("Oriental", "Occidental")

inx <- 2

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
    linewidth = ifelse(Type == "Masked", 1, 0.4)
  )

# Load color configuration for pouplation ----------------------------------
# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
ori_order_color <- setNames(color_uniq$Color, color_uniq$Population)

pops_in_data <- unique(counts_long$Population)
my_colors <- ori_order_color[pops_in_data]


# Combined plot
p1 <- ggplot(counts_long, aes(x = Cutoff, y = log10(Count),
                              group = interaction(Population, Type),
                              color = Population)) +
  geom_line(aes(linewidth = linewidth), linetype = "solid", alpha = 0.6) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color= "black",
            size = 0.4, alpha = 0.8)+
  geom_vline(xintercept = 0.005, linetype = "dashed", color= "black",
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

p1

# Faceted plot per population
p2 <- ggplot(counts_long, aes(x = Cutoff, y = log10(Count),
                              group = interaction(Population, Type),
                              color = Population)) +
  geom_line(aes(linewidth = linewidth), linetype = "solid", alpha = 0.6) +
  geom_point(size = 1, shape = 16) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color= "black",
             size = 0.4, alpha = 0.8)+
  geom_vline(xintercept = 0.005, linetype = "dashed", color= "black",
             size = 0.4, alpha = 0.8)+
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
  filename = paste0("../plot/s05.a2.TajimaD_GeneCount_Distribution_", folder[inx], ".png"),
  plot = final_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600
)

ggsave(
  filename = paste0("../plot/s05.a2.TajimaD_GeneCount_Distribution_", folder[inx], ".tif"),
  plot = final_plot,
  width = 180 / 25.4,
  height = 100 / 25.4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

