library(ggplot2)
library(ggpubr)
library(dplyr)
library(tune)
library("scatterplot3d") 
library("plotly")


rm(list=ls())
# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_path)


# Read PCA data
# Western 
prefix_Eu <- "../output/branch15.batch40.pear_Dec2023.Europe.remove_others/pear_Dec2023.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k.pca"
# Eastern
prefix_As <- "../output/branch15.batch41.Asia.remove_others/pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k.pca"

prefix <- prefix_As


pca_path <- paste0(prefix, ".eigenvec.input_PCA_plot.data")
data <- read.table(pca_path,sep ='\t', header = TRUE )


# Read groups and colors information data, and merge data
pop <- read.table("../input/population_group.tsv", sep = '\t', header = TRUE)
data <- left_join(data, select(pop, ID_vcf, Population, Color), by = "ID_vcf")

# Group mapping
group_mapping <- c("CE" = "Cultivated",
                   "WE" = "Wild",
                   "CW" = "Cultivated",
                   "WW" = "Wild")

# Population mapping
data <- data %>%
  mutate(Population = as.character(Population),
         Population = ifelse(substr(Population, nchar(Population) - 3, nchar(Population)) == "Admi", "Admixed", Population)
        )

# Change the value according to mappings
data <- data %>%
  mutate(Group = recode(Group, !!!group_mapping),
         # Group = ifelse(Population == "Admixed", Population, Group),
         Group = ifelse(Population == "betu_D3C1", "Rootstock", Group)
         )


# Using shiny to have an impression
plot_ly(data,x=data$PC1,y=data$PC3, z=data$PC2, color=data$Color)
# ?plot_ly


# Define the title
if (prefix == prefix_Eu) {
  pca_title <- "PCA for Western Pear"
} else if (prefix == prefix_As) {
  pca_title <- "PCA for Eastern Pear"
}

# Define the shapes 
if (prefix == prefix_Eu) {
  # shapes_index <- c(1,4,17)
  shapes_index <- c(4,17)   # Update on 2024-05-29 from the last line, as Admixed doesn't use specific shape
  pch_shapes_index <- shapes_index[c(2,1)]  #   This is used in plot legend
} else if (prefix == prefix_As) {
  #  shapes_index <- c(1,8,7,6)
  # pch_shapes_index <-  shapes_index[c(4,1,2,3)]  # For As dataset, The orders in legend changed a bit
  shapes_index <- c(4,7,17)  # Update on 2024-05-29 from the last line, as Admixed doesn't use specific shape
  pch_shapes_index <-  shapes_index[c(1,3,2)]  
}
shapes <- shapes_index[as.numeric(as.factor(data$Group))]


# FUN: Read the pve file
pve_path <- paste0(prefix,".pve")
read_pve <- function(pve_path, index) {
    pve <- read.csv(pve_path, header = FALSE)
    pve <- paste0(round(pve[index,], 1),"%")
    return(pve)}


# Plot 
#?scatterplot3d
main_plot <- scatterplot3d(data$PC1, data$PC2, data$PC3, 
              color=data$Color, 
              angle = 20,
              scale.y = 1,
              grid = T,
              box = F,
              pch = shapes,
              # main = pca_title,
              xlab = paste0("PC1(",read_pve(pve_path,1),")"),
              ylab = paste0("PC2(",read_pve(pve_path,2),")"),
              zlab = paste0("PC3(",read_pve(pve_path,3),")"),
              cex.lab=0.8,
              cex.axis=0.8,
              tick.marks = T,
              label.tick.marks = F
)


main_plot <- scatterplot3d(data$PC1, data$PC3, data$PC2, 
                           color=data$Color, 
                           angle = 10,
                           scale.y = 1,
                           grid = T,
                           box = F,
                           pch = shapes,
                           # main = pca_title,
                           xlab = paste0("PC1(",read_pve(pve_path,1),")"),
                           ylab = paste0("PC3(",read_pve(pve_path,3),")"),
                           zlab = paste0("PC2(",read_pve(pve_path,2),")"),
                           cex.lab=0.8,
                           cex.axis=0.8,
                           tick.marks = T,
                           label.tick.marks = F
)


df <- arrange(data,Population)
main_plot <- scatterplot3d(df$PC1, df$PC3, df$PC2, 
                           color=df$Color, 
                           angle = 10,
                           scale.y = 1,
                           grid = T,
                           box = F,
                           pch = shapes,
                           # main = pca_title,
                           xlab = paste0("PC1(",read_pve(pve_path,1),")"),
                           ylab = paste0("PC3(",read_pve(pve_path,3),")"),
                           zlab = paste0("PC2(",read_pve(pve_path,2),")"),
                           cex.lab=0.8,
                           cex.axis=0.8,
                           tick.marks = T,
                           label.tick.marks = F
)


# Legend for type (shape indexed)
main_plot <-  legend("topleft", 
       pch=pch_shapes_index,
       legend=unique(data$Group),   # Check if map
       cex = .5
       )
 #?legend

# Legend for species (color indexed)
# data$Color <- as.factor(data$Color)
# Check the color and Population mapping by visual inspection
legend_data <- unique(data[, c("Color", "Population")])
legend_data
uniq_color <- unique(data$Color)
uniq_color

if (prefix == prefix_Eu) {
  legend_sp <- c("Admixed",
                 "P. pyraster", 
                 "P. communis", 
                 "P. communis", 
                 "P. communis",
                 "P. communis",
                 "P. caucasica")
} else if (prefix == prefix_As) {
  legend_sp <- c("P. pyrifolia Sand pear",
                 "P. pyrifolia Sand pear",
                 "P. pyrifolia Sand pear",
                 "Admixed",
                 "P. pyrifolia Japanese pear",
                 "P. pashia",
                 "P. pyrifolia White pear",
                 "P. betulifolia",
                 "P. ussuriensis",
                 "P. betulifolia")
}

legend("topright", pch = 19,
       col=unique(data$Color),
       legend=legend_sp,
       cex= 0.5
       )




######################  END  ##########


# FUN: Join the pca eu_data with pop and color setting file
join_pca_pop <- function(data) {
  data <- left_join(data, select(pop, ID_vcf, Population, Color), by = "ID_vcf")
  return(data)
}


# FUN: Read the pve file
read_pve <- function(prefix, index) {
  pve <- read.csv(paste0(prefix,"pve"), header = FALSE)
  pve <-paste0(round(pve[index,], 1),"%")
  return(pve)
}


# FUN: Set Population names
process_population_data <- function(data) {
  data <- data %>%
    mutate(Population = as.character(Population),
           Population = ifelse(substr(Population, nchar(Population) - 3, nchar(Population)) == "Admi", "Admixed", Population))
  
  data <- data %>%
    group_by(Population) %>%
    mutate(Population_count = n()) %>%
    ungroup()
  
  data$Population <- ifelse(data$Population_count < 5, "Others", data$Population)
  
  return(data)
}


# FUN: Define the plot function
word_size <- theme(text = element_text(size = 8),  
                   axis.title = element_text(size = 7), 
                   axis.text = element_text(size = 6), 
                   plot.title = element_text(size = 7))

generate_3d_scatter_plot <- function(data, x_var,y_var,z_var, pve_x, pve_y, pve_z,shape_var) {
  plot <- scatterplot3d(x=x_var, y=y_var, z=z_var)
  
  # plot <- ggplot(data, aes(x = x_var, y = y_var, label = Group)) +
  #   geom_point(aes(shape = Group, color = Color), size = 1.5, alpha = 1) +
  #   xlab(paste0("PC1 (", pve_x, ")")) +
  #   ylab(paste0("PC2 (", pve_y, ")")) +
  #   theme_classic() +
  #   #theme_void()+
  #   coord_fixed(ratio = 1)+
  #   scale_shape_manual(values = shape_var) +
  #   scale_color_identity(name = "Population", 
  #                        labels = data$Population, 
  #                        breaks = data$Color, 
  #                        guide = "legend")+
  #   word_size +
  #   tune::coord_obs_pred()
  
  return(plot)
}





# FUN: Main function
main_plots <- function(prefix, data_suffix, shape){
  pve1 <- read_and_process_pve(prefix, 1)
  pve2 <- read_and_process_pve(prefix, 2)
  pve3 <- read_and_process_pve(prefix, 3)

  pop <- read.table("../input/population_group.tsv", sep = '\t', header = TRUE)
  
  data <- read.table(paste0(prefix,suffix),sep ='\t', header = TRUE)
  data <- join_pca_pop(data)
  data <- process_population_data(data)
  
  # Explain abbreviations in code
  group_mapping <- c("CE" = "Cultivated Eastern",
                     "WE" = "Wild Eastern",
                     "CW" = "Cultivated Western",
                     "WW" = "Wild Western")
  data <- data %>%
    mutate(Group = recode(Group, !!!group_mapping),
           Group = ifelse(Population == "Admixed", Population, Group))
  
  data$Group <- as.factor(data$Group)
  data$Population <- as.factor(data$Population)
  data$Color <- as.factor(data$Color)

  word_size <- theme(text = element_text(size = 8),  
                   axis.title = element_text(size = 7), 
                   axis.text = element_text(size = 6), 
                   plot.title = element_text(size = 7))
  
  p1 <- generate_3d-scatter_plot(data, data$PC1, data$PC2, data$PC3, pve1, pve2,pve3, shape)
  #p2 <- generate_scatter_plot(data, data$PC1, data$PC3, pve1, pve3, shape)
  #p3 <- generate_scatter_plot(data, data$PC2, data$PC3, pve2, pve3, shape)

  #plots <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = 'bottom')                 
}

prefix ='test.data'





################################################################################

#  Following codes are running without main function, just keep in the pocket.

################################################################################

# Read the pca data file
eu_prefix ='output/branch15.batch27.pear_Dec2023.Europe/pear_Dec2023.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k.pca.'
eu_data <- read.table(paste0(eu_prefix,"eigenvec.input_PCA_plot.K11.data"),
                   sep ='\t', header = TRUE)

as_prefix <- 'output/branch15.batch31.Asia.removeConflicts/pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k.pca.'
as_data <- read.table(paste0(as_prefix,"eigenvec.input_PCA_plot.K10.data"),
                   sep ='\t', header = TRUE)

# Read pop and color setting file
pop <- read.table("input/population_group.tsv", sep = '\t', header = TRUE)

# Read the pve file
eu_pve_1 <- read_and_process_pve(eu_prefix, 1)
eu_pve_2 <- read_and_process_pve(eu_prefix, 2)
eu_pve_3 <- read_and_process_pve(eu_prefix, 3)
eu_pve_4 <- read_and_process_pve(eu_prefix, 4)

as_pve_1 <- read_and_process_pve(as_prefix, 1)
as_pve_2 <- read_and_process_pve(as_prefix, 2)
as_pve_3 <- read_and_process_pve(as_prefix, 3)
as_pve_4 <- read_and_process_pve(as_prefix, 4)


# Join the pca eu_data with pop and color setting file
eu_data <- join_pca_pop(eu_data)
as_data <- join_pca_pop(as_data)

# Set Population names
eu_data <- process_population_data(eu_data)
as_data <- process_population_data(as_data)

# Explain abbreviations in code
# Define the mapping for Group column
group_mapping <- c("CE" = "Cultivated Eastern",
                   "WE" = "Wild Eastern",
                   "CW" = "Cultivated Western",
                   "WW" = "Wild Western")

# Apply the mapping to eu_data
eu_data <- eu_data %>%
  mutate(Group = recode(Group, !!!group_mapping),
         Group = ifelse(Population == "Admixed", Population, Group))

# Apply the mapping to as_data
as_data <- as_data %>%
  mutate(Group = recode(Group, !!!group_mapping),
         Group = ifelse(Population == "Admixed", Population, Group))


# Convert Group and Color columns to factors
eu_data$Group <- as.factor(eu_data$Group)
eu_data$Population <- as.factor(eu_data$Population)
eu_data$Color <- as.factor(eu_data$Color)

as_data$Group <- as.factor(as_data$Group)
as_data$Population <- as.factor(as_data$Population)
as_data$Color <- as.factor(as_data$Color)

eu_shapes <- c(1,4,17)
p_eu12 <- generate_scatter_plot(eu_data,eu_data$PC1, eu_data$PC2, eu_pve_1, eu_pve_2, eu_shapes)
p_eu13<- generate_scatter_plot(eu_data,eu_data$PC1, eu_data$PC3, eu_pve_1, eu_pve_3, eu_shapes)
p_eu23 <- generate_scatter_plot(eu_data,eu_data$PC2, eu_data$PC3, eu_pve_2, eu_pve_3, eu_shapes)


plot_eu <-ggarrange(p_eu12,p_eu13,p_eu23,ncol=3, nrow=1, common.legend = TRUE, legend = 'bottom')
plot_eu


as_shapes <- c(1, 16, 6)
p_as12 <- generate_scatter_plot(as_data,as_data$PC1, as_data$PC2, as_pve_1, as_pve_2, as_shapes)
p_as13<- generate_scatter_plot(as_data,as_data$PC1, as_data$PC3, as_pve_1, as_pve_3, as_shapes)
p_as23 <- generate_scatter_plot(as_data,as_data$PC2, as_data$PC3, as_pve_2, as_pve_3, as_shapes)

plot_as <-ggarrange(p_as12,p_as13,p_as23,ncol=3, nrow=1, common.legend = TRUE, legend = 'bottom')
plot_as







################################################################################

#  Following codes are unnecessary, just keep in the pocket.

################################################################################
p_eu12 <- ggplot(eu_data, aes(x = PC1, y = PC2,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC1 (", eu_pve_1, ")")) +
  ylab(paste0("PC2 (", eu_pve_2, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  # theme_linedraw()+
  # ylim(-0.2,0.3) +
  # xlim(-0.2,0.3) +
  scale_shape_manual(values = c(1, 4, 17)) +
  # scale_x_continuous(expand = c(0, 0), limits = c(-0.1, 0.35)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-0.15, 0.3)) +
  scale_color_identity(name = "Population", 
                       labels = eu_data$Population, 
                       breaks = eu_data$Color, 
                       guide = "legend")+
  word_size


p_eu12 <- p_eu12 + tune::coord_obs_pred()





p_eu13 <- ggplot(eu_data, aes(x = PC1, y = PC3,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC1 (", eu_pve_1, ")")) +
  ylab(paste0("PC3 (", eu_pve_3, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  scale_shape_manual(values = c(1, 4, 17)) +
  scale_color_identity(name = "Population", 
                       labels = eu_data$Population, 
                       breaks = eu_data$Color, 
                       guide = "legend")+
  word_size

p_eu13 <- p_eu13 + tune::coord_obs_pred()

p_eu23 <- ggplot(eu_data, aes(x = PC2, y = PC3,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC2 (", eu_pve_2, ")")) +
  ylab(paste0("PC3 (", eu_pve_3, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  scale_shape_manual(values = c(1, 4, 17)) +
  scale_color_identity(name = "Population", 
                       labels = eu_data$Population, 
                       breaks = eu_data$Color, 
                       guide = "legend")+
  word_size

p_eu23 <- p_eu23 + tune::coord_obs_pred()



plot_eu <- ggarrange(p_eu12,p_eu13,p_eu23,ncol=3, nrow=1, common.legend = TRUE, legend = 'bottom')
plot_eu 

ggsave("plot_eu.pdf", plot_eu, width = 85, height = 74.25, units = "mm", device = "pdf")
ggsave("plot_eu2.pdf", plot_eu, width = 210, height = 74.25, units = "mm", device = "pdf")


p_as12 <- ggplot(as_data, aes(x = PC1, y = PC2,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC1 (", as_pve_1, ")")) +
  ylab(paste0("PC2 (", as_pve_2, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  scale_shape_manual(values = c(1, 16, 6)) +
  scale_color_identity(name = "Population", 
                       labels = as_data$Population, 
                       breaks = as_data$Color, 
                       guide = "legend")+
  tune::coord_obs_pred()+
  word_size
p_as12


p_as13 <- ggplot(as_data, aes(x = PC1, y = PC3,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC1 (", as_pve_1, ")")) +
  ylab(paste0("PC3 (", as_pve_3, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  scale_shape_manual(values = c(1, 16, 6)) +
  scale_color_identity(name = "Population", 
                       labels = as_data$Population, 
                       breaks = as_data$Color, 
                       guide = "legend")+
  tune::coord_obs_pred()+
  word_size

p_as13

p_as23 <- ggplot(as_data, aes(x = PC2, y = PC3,  label = Group)) +
  geom_point(aes(shape = Group, color = Color), size = 0.8, alpha = 1) +
  xlab(paste0("PC2 (", as_pve_2, ")")) +
  ylab(paste0("PC3 (", as_pve_3, ")")) +
  theme_classic() +
  coord_fixed(ratio = 1)+
  scale_shape_manual(values = c(1, 16, 6)) +
  scale_color_identity(name = "Population", 
                       labels = as_data$Population, 
                       breaks = as_data$Color, 
                       guide = "legend")+
  tune::coord_obs_pred()+
  word_size
p_as23

plot_as <-ggarrange(p_as12,p_as13,p_as23,ncol=3, nrow=1, common.legend = TRUE, legend = 'bottom')
plot_as
ggsave("plot_as.pdf", plot_as, width = 85, height = 74.25, units = "mm", device = "pdf")
ggsave("plot_as2.pdf", plot_as, width = 210, height = 74.25, units = "mm", device = "pdf")



# ggarrange(p_as12, p_as13, p_as23, ncol = 3, common.legend = TRUE, legend = 'bottom')

legend_as <- get_legend(p_as12)
legend_eu <- get_legend(p_eu12)
legends <- ggarrange(legend_eu, legend_as, nrow=2)

rm_legend <- function(p){p + theme(legend.position = "none")}

plots_eu <- ggarrange(rm_legend(p_eu12), rm_legend(p_eu13), rm_legend(p_eu23), 
                      ncol = 3, nrow =1)

plots_as <- ggarrange(rm_legend(p_as12), rm_legend(p_as13), rm_legend(p_as23),
                      ncol = 3, nrow =1)
plots <- ggarrange(plots_eu, plots_as, ncol=1, nrow=2)

plots <- plots_eu <- ggarrange(rm_legend(p_eu12), rm_legend(p_eu13), rm_legend(p_eu23), 
                               rm_legend(p_as12), rm_legend(p_as13), rm_legend(p_as23),
                      ncol = 3, nrow =2)
                   
                
ggarrange(plots, legends, widths = c(0.8,0.2) )




          # common.legend = TRUE, legend = 'bottom')+
  # guides(color = guide_legend(title = "Population", override.aes = list(shape = NA)))

# 
# library(patchwork)
# 
# combined <- p_eu12+ p_eu13+p_eu23+p_as12+p_as13+p_as23+ 
#     plot_layout(guides = "collect") & theme(legend.position = "bottom")
# 
# combined

library(cowplot)

p1 <- p_eu12
p2 <- p_eu13
p3 <- p_eu23
p4 <- p_as12
p5 <- p_as13
p6 <- p_as23



# arrange the three plots in a single row
prow <- plot_grid( p1 + theme(legend.position="none"),
                   p4 + theme(legend.position="none"),
                   p5 + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("A", "B", "C"),
                   hjust = -1,
                   nrow = 1
)

# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
legend_b <- get_legend(p1 + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p



# library(ggrepel)
# # install.packages("ggforce")
# library(ggforce)
# library(adegenet)

# library(tidyverse)
# library(eu_data.table)
# 
# library(FactoMineR)



# Set an option for avoiding label overlaps
# options(ggrepel.max.overlaps = 300)

# Create the ggplot object
# p02 <- ggplot(eu_data, aes(x = PC1, y = PC2,color=Color, shape = Group,label = Group)) +
#   geom_point(size = 2, alpha = 0.2) +
#   xlab("PC1 (30.0%) explained") +
#   ylab("PC2 (7.8%) explained") +
#   theme_classic() +
#   scale_color_identity()+
#   coord_fixed()
# p02
# 
# 
# 
# p03 <- ggplot(eu_data, aes(x = PC3, y = PC4,color=Color, shape = Group,label = Group)) +
#   geom_point(size = 2, alpha = 0.5) +
#   xlab("PC3 (%) explained") +
#   ylab("PC4 (%) explained") +
#   theme_classic() +
#   scale_color_identity()+
#   coord_fixed()
# p03
# 

# Plot with label and segments
# Function to calculate central points for each group
# calculate_central_points <- function(group_eu_data) {
#   central_point <- colMeans(group_eu_data[, c("PC1", "PC2","PC3","PC4")])
#   return(eu_data.frame(
#     Group = unique(group_eu_data$Group), 
#     Central_PC1 = central_point[1],
#     Central_PC2 = central_point[2], 
#     Central_PC3 = central_point[3],
#     Central_PC4 = central_point[4],
#     Color_c = group_eu_data$Color[1]))
# }

# # Calculate central points for each group
# central_points <- do.call(rbind, by(eu_data, eu_data$Group, calculate_central_points))
# nonadmix_central_points <- central_points[central_points$Group != "Admixed", ]
# 
# # Merge central points back into the PCA eu_data
# eu_data <- merge(eu_data, central_points, by = "Group")

# # in following code, the label has been shown too many times, it would be change in p05
# p04 <- ggplot(eu_data, aes(x = PC1, y = PC2,color=Color, shape = Group,label = Group)) +
#   geom_point(size = 2, alpha = 0.5) +
#   xlab("PC1 (%) explained") +
#   ylab("PC2 (%) explained") +
#   theme_classic() +
#   scale_color_identity()+
#   coord_fixed()+
#   geom_segment(aes(x = Central_PC1, y = Central_PC2, xend = PC1, yend = PC2), alpha = 0.5) +
#   geom_label(aes(x = Central_PC1, y = Central_PC2, label = Group), alpha = 0.2)
# p04




# # Filter out rows where Species is "Admix"
# nonadmix_eu_data <- eu_data[eu_data$Group != "Admixed", ]






# plot
concat <- cbind.eu_data.frame(eu_data$Group,eu_data$PC1,eu_data$PC2)
concat$`eu_data$Group` <- as.factor(concat$`eu_data$Group`)
ellipse.coord = coord.ellipse(concat,bary=T)
ellipse.coord <- as.eu_data.frame(ellipse.coord) 
setDT(ellipse.coord)
setnames(ellipse.coord,c("res.eu_data.Group","res.eu_data.PC1","res.eu_data.PC2"), c("Group", "x", "y"))

p05 <- ggplot(eu_data, aes(x = PC1, y = PC2,color=Group, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC1 (", pve_1, ") explained")) +
  ylab(paste0("PC2 (", pve_2, ") explained")) +
  theme_classic() +
  # scale_color_identity()+
  coord_fixed( ratio = 1)+
  geom_label(eu_data = nonadmix_central_points, aes(x = Central_PC1, y = Central_PC2, label =Group, color=Group), alpha = 0.5) +

geom_path(eu_data=ellipse.coord, aes(x=x, y=y, color = "black"),show.legend = FALSE)
print(p05)


concat <- cbind.eu_data.frame(eu_data$Color,eu_data$PC1,eu_data$PC2)
concat$`eu_data$Color` <- as.factor(concat$`eu_data$Color`)
ellipse.coord = coord.ellipse(concat,bary=T)
ellipse.coord <- as.eu_data.frame(ellipse.coord) 
setDT(ellipse.coord)
setnames(ellipse.coord,c("res.eu_data.Color","res.eu_data.PC1","res.eu_data.PC2"), c("Color", "x", "y"))

p051 <- ggplot(eu_data, aes(x = PC1, y = PC2,color=Color)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC1 (", pve_1, ") explained")) +
  ylab(paste0("PC2 (", pve_2, ") explained")) +
  theme_classic() +
  scale_color_identity()+
  coord_fixed( ratio = 1)+
  # geom_label(eu_data = nonadmix_central_points, aes(x = Central_PC1, y = Central_PC2, label =Group, color=Group), alpha = 0.5) +
  geom_path(eu_data=ellipse.coord, aes(x=x, y=y, color = 'black'),show.legend = FALSE)
print(p051)

# p05_1 <- p05 + geom_segment(eu_data=nonadmix_eu_data, aes(x = Central_PC1, y = Central_PC2, xend = PC1, yend = PC2), alpha = 0.5)



# Plot PC1 and PC3
p06 <- ggplot(eu_data, aes(x = PC1, y = PC3,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1 (30.0%) explained") +
  ylab("PC3 (4.8%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed( ratio = 1)+
  geom_label(eu_data = nonadmix_central_points, aes(x = Central_PC1, y = Central_PC3, label =Group, color=Color_c), alpha = 0.5)


# Plot PC 2 and PC3

# Plot with lable and segments
# 
p07 <- ggplot(eu_data, aes(x = PC2, y = PC3,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC2 (7.8%) explained") +
  ylab("PC3 (4.8%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()+
  geom_label(eu_data=nonadmix_central_points, aes(x = Central_PC2, y = Central_PC3, label =Group, color=Color_c), alpha = 0.5)

# p07_1 <- p07 + geom_segment(eu_data=nonadmix_eu_data, aes(x = Central_PC2, y = Central_PC3, xend = PC2, yend = PC3), alpha = 0.5)+

# print(p07)



# Plot PC 3 and PC4
p08 <- ggplot(eu_data, aes(x = PC3, y = PC4,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC3 (4.8%) explained") +
  ylab("PC4 (3.1%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed( ratio = 1) +
  geom_label(eu_data=nonadmix_central_points, aes(x = Central_PC3, y = Central_PC4, label =Group, color=Color_c), alpha = 0.5)

# p08_1 <- p08 + geom_segment(eu_data=nonadmix_eu_data, aes(x = Central_PC3, y = Central_PC4, xend = PC3, yend = PC4), alpha = 0.5)

ggarrange(p05, p06, p07, ncol=3, common.legend = TRUE,legend = 'bottom')

# ?ggarrange







pca_result <- prcomp(eu_data[c(5,24)], scale = TRUE)
summary(pca_result)$importance

scaled_eu_data <- scale(eu_data[c(5:24)])

eu_data <- cbind(eu_data, scaled_eu_data )

p07 <- ggplot(eu_data, aes(x = PC3.1, y = PC4.1,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC3 (%) explained") +
  ylab("PC4 (%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()
p07




# Following code is from mapmixture, the function of calculating centroid is the same with before

# # Calculate centroid (mean average) positions for each group
# centroid_df <- stats::aggregate(
#   cbind(PC1, PC2) ~ Group,
#   FUN = mean,
#   eu_data = eu_data
# )
# 
# # Edit column names of centroid eu_data frame
# colnames(centroid_df)[2:3] <- c("cen1", "cen2")
# 
# # Add centroid coordinates to main eu_data frame
# df <- dplyr::left_join(eu_data, centroid_df, by = "Group")

