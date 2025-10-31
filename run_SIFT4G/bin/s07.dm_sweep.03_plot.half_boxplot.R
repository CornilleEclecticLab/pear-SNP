# Refered Xilong Chen's script

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(multcompView)

set.seed(42) # for jitter

# Set working directory to the location of the script and clean up
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
  rm(list=ls())
}


# Function save plot, w-max = 180mm, h-max = 170mm
save_myplot <- function(p, prefix, w = 180, h = 100, dpi = 600) {
  wi <- w / 25.4
  hi <- h / 25.4
  
  ggsave(paste0(prefix, ".pdf"), p,
         width = wi, height = hi, units = "in")
  
  ggsave(paste0(prefix, ".tif"), p,
         width = wi, height = hi, units = "in",
         dpi = dpi, device = "tiff", compression = "lzw")
  
  ggsave(paste0(prefix, ".png"), p,
         width = wi, height = hi, units = "in", dpi = dpi)
}


# Function transfer triangle df to square df 
tri.to.squ <- function(tri) {
  rows <- rownames(tri)
  cols <- colnames(tri)
  all_groups <- union(rows, cols)
  
  mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                dimnames = list(all_groups, all_groups))
  
  for (i in rows) {
    for (j in cols) {
      val <- tri[i, j]
      if (!is.na(val)) {
        mat[i, j] <- val
        mat[j, i] <- val
      }
    }
  }
  
  return(mat)
}



# 08-22
### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)

# Define population order for Occidental/Eastern
# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "pyra", "cauc",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]
wild_pops     <- c("pyra", "cauc", "pash", "ussu", "betu")

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels),
         region_type = factor(region_type, levels = c("all_masked", "control", "selection")),
         crop_type = if_else(pop %in% wild_pops, "Wild", "Cultivar"))

df_sweep$de_per_100k_cds <- df_sweep$de_per_cds * 100000


# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
pop_color <- setNames(color_uniq$Color, color_uniq$Population)

# pops_in_data <- unique(counts_long$Population)
# my_colors <- ori_order_color[pops_in_data]

