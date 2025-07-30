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

input_dir_real <- file.path(base_dirs[inx], "output/s02.split_pixy_by_pop")

# List all reports in the directory
file_list_real <- list.files(path = input_dir_real, pattern = "^Pixy\\..*\\.1k_tajima_d\\.txt$", full.names = TRUE)

# Function to read and clean TajimaD real results
read_real_values <- function(file) {
  print(paste0("Processing ", file))  
  file_name <- basename(file)
  pop_chr <- str_split(file_name, "\\.", simplify = TRUE)
  population <- pop_chr[2]
  chromosome <- pop_chr[3]
  
  # Read and preprocess the file
  file_lines <- readLines(file)
  file_lines <- file_lines[!(file_lines == "" | startsWith(file_lines, "//"))]  # Drop empty lines and lines starting with //
  
  # header:
  # pop chromosome window_pos_1 window_pos_2 tajima_d no_sites raw_pi raw_watterson_theta tajima_d_stdev
  df <- read.table(file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

  df <- df %>%
    drop_na(tajima_d) %>%          # rm NA
    filter(is.finite(tajima_d))    # rm nan/Inf

  df$Population <- population
  df$Chromosome <- chromosome
  df$Dataset <- "Real"

  df <- df %>% 
    mutate(tajima_d = tajima_d)
  return(df)
}

# Read all real files and combine data
data_real <- bind_rows(lapply(file_list_real, read_real_values))


# Load color configuration for pouplation ----------------------------------
# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
ori_order_color <- setNames(color_uniq$Color, color_uniq$Population)

# Purge color
populations_in_data <- unique(data_real$Population)
my_colors <- c(ori_order_color[names(ori_order_color) %in% populations_in_data])




# Calculation and Plot -----------------------------------------------------

# Add log10(TajimaD) column
# No, beacuse of the values < 0

# data_real <- data_real %>%
#   mutate(LogDTajimaD= log10(TajimaD)) %>%
#   filter(is.finite(LogTajimaD))


# Compute quantiles and prepare for plotting & labeling
q01  <- data_real %>% group_by(Population) %>% summarise(cut = quantile(tajima_d, 0.01),  .groups = "drop") %>% mutate(tail = "1%")
q005 <- data_real %>% group_by(Population) %>% summarise(cut = quantile(tajima_d, 0.005), .groups = "drop") %>% mutate(tail = "0.5%")
q001 <- data_real %>% group_by(Population) %>% summarise(cut = quantile(tajima_d, 0.001), .groups = "drop") %>% mutate(tail = "0.1%")


qall <- bind_rows(q01, q005, q001) %>%
  # assign line properties and y-position for labels
  mutate(
    linetype  = ifelse(tail == "0.5%", "solid", "dashed"),
    linewidth = ifelse(tail == "0.1%", 0.3, 0.2),
    y_pos     = runif(n(), 0.2, 0.8)
  )


# Compute mean, sd, cutoffs, zscore and p-value per population
stats <- data_real %>% 
  group_by(Population) %>% 
  summarise(
    mean_td = mean(tajima_d),
    sd_td   = sd(tajima_d),
    cut_01  = quantile(tajima_d, 0.01),
    cut_005 = quantile(tajima_d, 0.005),
    cut_001 = quantile(tajima_d, 0.001),
    .groups = "drop"
  ) %>% 
  mutate(
    z_01  = (cut_01  - mean_td) / sd_td,
    p_01  = pnorm(z_01),          # bottom/left tail
    z_005 = (cut_005 - mean_td) / sd_td,
    p_005 = pnorm(z_005),
    z_001 = (cut_001 - mean_td) / sd_td,
    p_001 = pnorm(z_001)
  )

stats_long <- stats %>% 
  select(Population, z_01, z_005, z_001) %>% 
  pivot_longer(cols = starts_with("z_"),
               names_to = "tail", values_to = "zscore") %>% 
  mutate(
    tail = recode(tail,
                  z_01  = "1%",
                  z_005 = "0.5%",
                  z_001 = "0.1%")
  )

qall <- left_join(qall, stats_long, by = c("Population", "tail"))



qd01  <- data_real %>%
  group_by(Population) %>%
  summarise(cut = quantile(tajima_d, 0.01), .groups="drop") %>%
  mutate(tail="99%")
qd001 <- data_real %>%
  group_by(Population) %>%
  summarise(cut = quantile(tajima_d, 0.001), .groups="drop") %>%
  mutate(tail="99.9%")


# Compute mean tajima_d per population
mean_df <- data_real %>% 
  group_by(Population) %>% 
  summarise(mean_tajima_d = mean(tajima_d), .groups = "drop") %>% 
  mutate(label_y = runif(n(), 0.2, 0.8),
         label_text = sprintf("Mean: %.2f", mean_tajima_d))



# # Function to count genes
# count_genes <- function(path) {
#   if (!file.exists(path)) return(NA_integer_)
#   length(readLines(path))
# }

# # Build outlier counts table
# populations <- unique(data_real$Population)
# counts_df <- bind_rows(lapply(populations, function(pop) {
#   # Raw outliers (s04) and masked outliers (s05)
#   f99  <- file.path(base_dirs[inx], "output", "s04.catch_outliers",
#                     sprintf("RAiSD.0_01outliers.%s.merged.bed.genes.txt", pop))
#   f999 <- file.path(base_dirs[inx], "output", "s04.catch_outliers",
#                     sprintf("RAiSD.0_001outliers.%s.merged.bed.genes.txt", pop))
#   f9999 <- file.path(base_dirs[inx], "output", "s04.catch_outliers",
#                     sprintf("RAiSD.0_0001outliers.%s.merged.bed.genes.txt", pop))
  
#   m99  <- file.path(base_dirs[inx], "output", "s05.mask_outliers",
#                     sprintf("RAiSD.0_01outliers.%s.genes.txt", pop))
#   m999 <- file.path(base_dirs[inx], "output", "s05.mask_outliers",
#                     sprintf("RAiSD.0_001outliers.%s.genes.txt", pop))
#   m9999 <- file.path(base_dirs[inx], "output", "s05.mask_outliers",
#                     sprintf("RAiSD.0_0001outliers.%s.genes.txt", pop))
  
#   tibble(
#     Population = pop,
#     tail       = c("99%", "99.9%", "99.99%"),
#     cut        = c(q99$cut[q99$Population==pop], q999$cut[q999$Population==pop], q9999$cut[q9999$Population==pop]),
#     n_raw      = c(count_genes(f99),  count_genes(f999), count_genes(f9999)),
#     n_mask     = c(count_genes(m99),  count_genes(m999), count_genes(m9999))
#   )
# }))

# # Prepare count labels
# counts_df <- counts_df %>%
#   left_join(
#     qall %>% select(Population, tail, y_pos),
#     by = c("Population", "tail")
#   ) %>%
#   mutate(
#     label2 = paste0("#Raw", n_raw, " Masked", n_mask))

data_r <- data_real
# data_real <- data_r %>% filter(Population=='pash')


# Plot 1: combined
p1 <- ggplot(data_real, aes(x=tajima_d, color=Population)) +
  geom_density(linewidth=0.3, alpha=0.5) +
  geom_rug(sides="b", alpha=0.3) +
  geom_vline(data=qall, aes(xintercept=cut, linetype=linetype, linewidth=linewidth, color=Population),
             show.legend=FALSE) +
  
  # Labels for cutoffs
  geom_text(
    data=qall,
    aes(x=cut, y=y_pos, label=paste0(tail,": ", sprintf("%.2f", cut), " Z=", sprintf("%.2f", zscore) ), color=Population),
    angle=90, vjust=-0.4, size=1, family="Arial", show.legend=FALSE
  ) +
  # Outlier gene count labels
  # geom_text(
  #   data=counts_df,
  #   aes(x=cut, y=y_pos, label=label2, color=Population),
  #   angle=90, vjust=1.2, size=1, family="Arial", show.legend=FALSE
  # ) +
  geom_vline(data = mean_df, aes(xintercept = mean_tajima_d, color = Population),
               linetype = "dotdash", linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = mean_df, aes(x = mean_tajima_d, y = label_y, label = label_text, color = Population),
            angle = 90, vjust = -0.4, size = 1, family = "Arial", show.legend = FALSE) +
  scale_color_manual(values=my_colors) +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  labs(title=paste0("Combined TajimaD Density with Cutoffs(", folder[inx], ")"),
       x="TajimaD", y="Density") +
  theme_minimal(base_size=6) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    axis.text = element_text(size = 5, family = "Arial"),
    legend.title = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 5, family = "Arial")
  )# +
  #coord_cartesian(xlim = c(-5, 5))  

# Plot 2: faceted by population
p2 <- ggplot(data_real, aes(x = tajima_d, color = Population)) +
  geom_density(linewidth=0.3, alpha = 0.5) +
  geom_rug(sides = "b", alpha = 0.3) +
  geom_vline(data = qall, aes(xintercept = cut, linetype = linetype, linewidth = linewidth, color = Population),
             show.legend = FALSE) +
  geom_text(
    data=qall,
    aes(x=cut, y=y_pos, label=paste0(tail,": ", sprintf("%.2f", cut), " Z=", sprintf("%.2f", zscore)), color=Population),
    angle=90, vjust=-0.4, size=1, family="Arial", show.legend=FALSE
  ) +
  # geom_text(
  #   data=counts_df,
  #   aes(x=cut, y=y_pos, label=label2, color=Population),
  #   angle=90, vjust=1.2, size=1, family="Arial", show.legend=FALSE
  # ) +
  geom_vline(data = mean_df, aes(xintercept = mean_tajima_d, color = Population),
               linetype = "dotdash", linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = mean_df, aes(x = mean_tajima_d, y = label_y, label = label_text, color = Population),
            angle = 90, vjust = -0.4, size = 1, family = "Arial", show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  facet_wrap(~Population, nrow = 4, ncol = 4) +
  labs(title = "Faceted TajimaD Density with Cutoffs per Population", 
       x = "TajimaD", y = "Density") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    axis.text = element_text(size = 5, family = "Arial"),
    strip.text = element_text(size = 6, family = "Arial"),
    legend.position = "none"
  ) #+
  #coord_cartesian(xlim = c(-5, 5))  

# Combine plots
final_plot <- (p1 / p2) + plot_layout(widths = c(1, 1), guides = "collect")
print(final_plot)


# Set plot size
width_in <- 180 / 25.4  # 180 mm to ~7.09 inches)
height_in <- 100 / 25.4  # 100 mm to ~3.94 inches)
dpi <- 600
pixel_width <- width_in * dpi
pixel_height <- height_in * dpi


# 
ggsave(
  filename = paste0("preview_TajimaD_10k_Destribution_unlimit_percent_", folder[inx], ".png"),
  plot = final_plot,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = dpi,
  device = "png"
)

# Export as TIFF
ggsave(
  filename = paste0("preview_TajimaD_10k_Destribution_unlimit_percent_", folder[inx], ".tif"),
  plot = final_plot,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = dpi,
  device = "tiff",
  compression = "lzw"
)





