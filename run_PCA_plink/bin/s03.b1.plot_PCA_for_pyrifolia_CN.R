# This is a clean version of the PCA plot script, which is used to generate 2D PCA plots for the data.
# When writing this script, referenced previous s03.plot_PCA.R script and s04.plot_3D_PCA.R, and removed all the unnecessary code.
# Date: 2024-07-31 16:04:27


# Load libraries
library(ggplot2)
library(dplyr)
library(rstudioapi)


# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
show(script_path)
setwd(script_path)


# Clean the environment
rm(list=ls())


# Set data path and load data
# pyri_CN
prefix_pyri_CN <- "../output/branch15.batch37.Cultivar_pyri_CN/pear_Dec2023.Cultivar_pyri_CN.Combine_chr.geno20_maf005.anno.syno.thin8k.pca"
prefix <- prefix_pyri_CN

pca_path <- paste0(prefix, ".eigenvec.input_PCA_plot.data")
data <- read.table(pca_path,sep ='\t', header = TRUE )


# Read groups and colors information data, and merge data
pop <- read.table("../input/population_group.for_Eastern_and_Western.2024-07-29.tsv", sep = '\t', header = TRUE)
data <- left_join(data, select(pop, ID_vcf, Population, Color), by = "ID_vcf")


# Group mapping
group_mapping <- c("CE" = "Cultivated",
                   "WE" = "Wild",
                   "CW" = "Cultivated",
                   "WW" = "Wild",
                   "RE" = "Rootstock")


# Population mapping
data <- data %>%
  mutate(Population = as.character(Population),
         Population = ifelse(substr(Population, nchar(Population) - 3, nchar(Population)) == "Admi", "Admixed", Population)
  )


# Change the value according to mappings
data <- data %>%
  mutate(Group = recode(Group, !!!group_mapping)    
  )


# FUN: Read the pve file
pve_path <- paste0(prefix,".pve")
read_pve <- function(pve_path, index) {
  pve <- read.csv(pve_path, header = FALSE)
  pve <- paste0(round(pve[index,], 1),"%")
  return(pve)
}


# Custom shapes for each group
shape_mapping <- c("Cultivated" = 4,  # 
                   "Wild" = 17,       # Triangle
                   "Rootstock" = 7)   # 


# Define the word size
word_size <- theme(text = element_text(size = 8),  
                   axis.title = element_text(size = 7), 
                   axis.text = element_text(size = 6), 
                   plot.title = element_text(size = 7))


# Plot PC1 vs PC2
ellipse_df <- subset(data, Population != "Admixed" & Population != "Others")

plot <- ggplot(data, aes(x = PC1, y = PC2, label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC1 (", read_pve(pve_path, 1), ")")) +
  ylab(paste0("PC2 (", read_pve(pve_path, 2), ")")) +
  theme_classic() +
  coord_fixed(ratio = 1) +
  scale_color_identity(name = "Population",
                       labels = unique(data$Population),
                       breaks = unique(data$Color),
                       guide = "legend"
  ) +
  scale_shape_manual(name = "Group", values = shape_mapping) +
  # Add ellipses
  stat_ellipse(data = ellipse_df, aes(group = Population, color = Color), level = 0.95, alpha = 1, size = 0.2) +
  word_size

# Print the plot
print(plot)
ggsave(paste0(prefix, ".PCA.PC1_vs_PC2.pdf"), width = 4, height = 4)