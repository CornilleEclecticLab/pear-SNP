library(ggplot2)
library(gplots)
library(RColorBrewer)


rm(list=ls())
# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_path)

data <- read.table("../output/s03.calculate_Pi_Dxy_across_chromosome/Fst_Pixy_in_Upper_and_Dxy_Pixy_in_Lower_matrix.tsv",
                   header=T, 
                   sep='\t', 
                   row.names=1, 
                   na.strings = "-")

data <- as.matrix(data)


# # Method 1 heatmap 
# heatmap(data, Colv = "Rowv", Rowv = NA)
# 
# my_group <- c(rep(1, 4), rep(2, 8))
# colSide <- brewer.pal(9, "Set1")[my_group]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data , RowSideColors=colSide , Colv = "Rowv", Rowv = NA)
# 
# # Method 2 heatmap.2
# library(gplots)
# heatmap.2(data)


# Method 3 pheatmap
display.brewer.all(colorblindFriendly = TRUE)
library("pheatmap")

# CRAN version
# install.packages("cartography")
library(cartography)

# logged_data <- log10(data)
# ?log

# upper_palette <- colorRampPalette(rev(brewer.pal(n = 4, name = "Greens")))(100)
# lower_palette <- colorRampPalette(rev(brewer.pal(n = 3, name = "BuPu")))(100)

upper_palette <- colorRampPalette(carto.pal(n1 = 8, pal1 = "kaki.pal"))(100)
lower_palette <- colorRampPalette(c( "white", "#555555"))(100)



upper_data <- data
upper_data[lower.tri(upper_data)] <- NA

pheatmap(upper_data, 
         color = upper_palette,
         cluster_rows = F, 
         cluster_cols = F,
         na_col = "White",
         annotation_colors = anno_colors,
         annotation_legend = F,
         number_format = "%0.4f",
         display_numbers = T,
         number_color = 'white'
)


lower_data <- data
lower_data[upper.tri(lower_data)] <- NA

pheatmap(lower_data, 
         color = lower_palette,
         cluster_rows = F, 
         cluster_cols = F,
         na_col = "White",
         annotation_colors = anno_colors,
         annotation_legend = F,
         number_format = "%0.4f",
         display_numbers = T,
         number_color = 'black'
)
 dev.size()


?pheatmap


