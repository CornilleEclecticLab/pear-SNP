suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(data.table)
  # library(parallel)
})

rm(list=ls())

# Set working directory to the location of the script (if running interactively)
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}


binary_cache_file <- "combined_data_binary.rds"
load(binary_cache_file)
print("Data loaded from binary file")

combined_data <- combined_data[combined_data$Score >= 1, ]

combined_data <- combined_data %>% 
  filter(!Population %in%  c("pash", "Sand_CN"))

head(combined_data)

############## TEST ################
# set.seed(123)
# sampled_data <- combined_data[sample(nrow(combined_data), 10000), ]
# combined_data <- sampled_data
####################################

#### Calculate boundry and midpoints of chromosome for labeling ####
chr_label <- combined_data %>%
  group_by(Chr) %>%
  summarize(
    end = max(Cum_pos)  # The maximum cumulative position for this chromosome
  ) %>%
  mutate(
    start = lag(end, default = 1),  # Start is the last max Cum_pos (default 1 for the first chromosome)
    midpoints = (start + end) / 2   # Calculate the midpoint for labeling
  )


# Load Color Configuration ------------------------------------------------------

# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>%
  distinct(Population, Color)
ori_order_color <- setNames(color_uniq$Color, color_uniq$Population)

# Purge color
populations_in_data <- unique(combined_data$Population)
my_colors <- c(ori_order_color[names(ori_order_color) %in% populations_in_data])
print("Loading color config")
my_colors

my_colors_with_grey <- c(my_colors, "Below_Cutoff" = "#AAAAAA")



# # Manhattan with dual y-axes ----------------------------------------------
# #### My theme
# my_theme <- theme_minimal() +
#   theme(
#     plot.title = element_blank(), # element_text(size = 8, face = "bold"),
#     plot.subtitle = element_blank(), # element_text(size = 8),
#     axis.title.x = element_blank(),# element_text(size = 6),
#     axis.title.y = element_text(size = 6),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#     axis.text.y = element_text(size = 5),
#     legend.title = element_text(size = 6),
#     legend.text = element_text(size = 6),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.ticks.x = element_line(),
#     axis.ticks.length = unit(.5, "mm"),
#     legend.position = "top",
#     panel.border = element_blank()
#   )
# 
# 
# #### Plot function
# plot_manhattan <- function(data, cutoff_data, title = "") {
#   ggplot(data) +
#     geom_point(
#       aes(x = Cum_pos, y = Score, color = Color_group, shape = Method),
#       alpha = .8, size = .5
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = raisd, color = Population),
#       linetype = "dashed", alpha = 1, linewidth = 0.2
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = omegaplus, color = Population),
#       linetype = "solid", alpha = 1, linewidth = 0.2
#     ) +
#     scale_x_continuous(
#       # name = "Chromosome", 
#       breaks = chr_label$midpoints,
#       labels = chr_label$Chr,
#       expand = c(0.01, 0.01)
#     ) +
#     scale_y_continuous(
#       name = "Mu",
#       expand = c(0.01, 0.01),
#       sec.axis = dup_axis(
#         # transform = ~ . / max(data$Score) * 12,  # Example transformation function
#         name = "Omega")
#     ) +
#     scale_color_manual(
#       values = my_colors_with_grey,
#       breaks = names(my_colors),
#       guide = guide_legend(override.aes = list(shape = 15, size = 1))
#     ) +
#     geom_vline(
#       xintercept = c(chr_label$start[1], chr_label$end), 
#       color = "grey50", linetype = "solid", alpha = 0.5
#     ) +
#     my_theme +
#     ggtitle(title)
# }
# 
# #### Plot p1, p2
# ## subset data
# cultivar_list <- c("comm_Dessert", "comm_Perry")
# non_cultivar_list <- c("cauc", "pyra")
# 
# if (!all(cultivar_list %in% unique(combined_data$Population))) {
#   cultivar_list <- c("White", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP")
#   non_cultivar_list <- c("ussu", "betu")
# }
# 
# sub_combined_data <- combined_data %>% 
#   filter(Population %in% cultivar_list)
# 
# sub_cutoff <- cutoff_df %>% 
#   filter(Population %in% cultivar_list)
# 
# non_cultivar_data <- combined_data %>% 
#   filter(Population %in% non_cultivar_list) %>%
#   mutate(Color_group = ifelse(Score > Cutoff, Population, "Below_Cutoff"))
# 
# non_cultivar_cutoff <- cutoff_df  %>% 
#   filter(Population %in% non_cultivar_list)
# 
# ## Plot
# print("Start plotting...")
# p1 <- plot_manhattan(sub_combined_data, sub_cutoff, "Cultivar Populations")
# p2 <- plot_manhattan(non_cultivar_data, non_cultivar_cutoff, "Non-Cultivar Populations")
# 
# 
# combined_plot <- plot_grid(p1, p2, ncol = 1)
# 
# combined_plot
# 
# 
# # Save
# ggsave("s07.plot_Manhattan.style2/combined_plot.pdf", plot = combined_plot, width = 180, height = 80, units = "mm", dpi = 300)
# ggsave("s07.plot_Manhattan.style2/combined_plot.png", plot = combined_plot, width = 180, height = 80, units = "mm", dpi = 300)
# ggsave("s07.plot_Manhattan.style2/combined_plot.A4.pdf", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 300)
# ggsave("s07.plot_Manhattan.style2/combined_plot.A4.png", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 300)
# 
# 
# 
# # Manhattan with dual y-axes facet_wrap by Pop  ----------------------------------------------
# 
# ## My theme
# my_theme <- theme_minimal() +
#   theme(
#     plot.title = element_blank(), # element_text(size = 8, face = "bold"),
#     plot.subtitle = element_blank(), # element_text(size = 8),
#     axis.title = element_text(size = 6),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#     axis.text.y = element_text(size = 5),
#     legend.title = element_text(size = 6),
#     legend.text = element_text(size = 6),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.ticks.x = element_line(),
#     axis.ticks.length = unit(.5, "mm"),
#     legend.position = "top",
#     panel.border = element_blank()
#   )
# 
# 
# ## Plot function
# plot_manhattan <- function(data, cutoff_data, title = "") {
#   ggplot(data) +
#     geom_point(
#       aes(x = Cum_pos, y = Score, color = Color_group, shape = Method),
#       alpha = .8, size = .5
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = raisd, color = Population),
#       linetype = "dashed", alpha = 1, linewidth = 0.2
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = omegaplus, color = Population),
#       linetype = "solid", alpha = 1, linewidth = 0.2
#     ) +
#     scale_x_continuous(
#       # name = "Chromosome", 
#       breaks = chr_label$midpoints,
#       labels = chr_label$Chr,
#       expand = c(0.01, 0.01)
#     ) +
#     scale_y_continuous(
#       name = "Mu",
#       expand = c(0.01, 0.01),
#       sec.axis = dup_axis(
#         # transform = ~ . / max(data$Score) * 12,  # Example transformation function
#         name = "Omega")
#     ) +
#     scale_color_manual(
#       values = my_colors_with_grey,
#       breaks = names(my_colors),
#       guide = guide_legend(override.aes = list(shape = 15, size = 1))
#     ) +
#     geom_vline(
#       xintercept = c(chr_label$start[1], chr_label$end), 
#       color = "grey50", linetype = "solid", alpha = 0.5
#     ) +
#     # add facet wrap
#     facet_wrap(~ Population, scales = "free_y") +
#     my_theme +
#     ggtitle(title)
# }
# 
# ## Plot p1, p2
# ## subset data
# cultivar_list <- c("comm_Dessert", "comm_Perry")
# non_cultivar_list <- c("cauc", "pyra")
# 
# if (!all(cultivar_list %in% unique(combined_data$Population))) {
#   cultivar_list <- c("White", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP")
#   non_cultivar_list <- c("ussu", "betu")
# }
# 
# sub_combined_data <- combined_data %>% 
#   filter(Population %in% cultivar_list)
# 
# sub_cutoff <- cutoff_df %>% 
#   filter(Population %in% cultivar_list)
# 
# non_cultivar_data <- combined_data %>% 
#   filter(Population %in% non_cultivar_list) %>%
#   mutate(Color_group = ifelse(Score > Cutoff, Population, "Below_Cutoff"))
# 
# non_cultivar_cutoff <- cutoff_df  %>% 
#   filter(Population %in% non_cultivar_list)
# 
# ## Plot
# print("Start plotting...")
# p1 <- plot_manhattan(sub_combined_data, sub_cutoff, "Cultivar Populations")
# p2 <- plot_manhattan(non_cultivar_data, non_cultivar_cutoff, "Non-Cultivar Populations")
# 
# 
# combined_plot <- plot_grid(p1, p2, ncol = 1)
# 
# combined_plot
# 
# 
# 
# # Manhattan with dual y-axes facet_wrap by Pop Chr 1 only without grid_plot----------------------------------------------
# 
# #### My theme 
# my_theme <- theme_minimal() +
#   theme(
#     plot.title = element_blank(), # element_text(size = 8, face = "bold"),
#     plot.subtitle = element_blank(), # element_text(size = 8),
#     axis.title = element_text(size = 6),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#     axis.text.y = element_text(size = 5),
#     legend.title = element_text(size = 6),
#     legend.text = element_text(size = 6),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.ticks.x = element_line(),
#     axis.ticks.length = unit(.5, "mm"),
#     legend.position = "top",
#     panel.border = element_blank()
#   )
# 
# 
# #### Plot function
# plot_manhattan <- function(data, cutoff_data, title = "") {
#   ggplot(data) +
#     geom_point(
#       aes(x = Pos, y = Score, color = Color_group, shape = Method),
#       alpha = .8, size = .5
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = raisd, color = Population),
#       linetype = "dashed", alpha = 1, linewidth = 0.2
#     ) +
#     geom_hline(
#       data = cutoff_data,
#       aes(yintercept = omegaplus, color = Population),
#       linetype = "solid", alpha = 1, linewidth = 0.2
#     ) +
#     scale_x_continuous(
#       breaks = chr_label$midpoints[chr_label$Chr == "Ch1"],
#       labels = chr_label$Chr[chr_label$Chr == "Ch1"],
#       expand = c(0.01, 0.01)
#     ) +
#     scale_y_continuous(
#       name = "Mu",
#       expand = c(0.01, 0.01),
#       sec.axis = dup_axis(
#         name = "Omega"
#       )
#     ) +
#     scale_color_manual(
#       values = my_colors_with_grey,
#       breaks = names(my_colors),
#       guide = guide_legend(override.aes = list(shape = 15, size = 1))
#     ) +
#     geom_vline(
#       xintercept = c(chr_label$start[chr_label$Chr == "Ch1"][1], chr_label$end[chr_label$Chr == "Ch1"]), 
#       color = "grey50", linetype = "solid", alpha = 0.5
#     ) +
#     facet_wrap(~ Population, scales = "free_y") +
#     my_theme +
#     ggtitle(title)
# }
# 
# #### Plotting for Ch1 only
# ## subset data for Ch1
# combined_data_ch1 <- combined_data %>% 
#   filter(Chr == "Chr1")
# 
# cutoff_df_ch1 <- cutoff_df %>% 
#   filter(Population %in% unique(combined_data_ch1$Population))
# 
# ## Plot
# print("Start plotting...")
# p <- plot_manhattan(combined_data_ch1, cutoff_df_ch1, "Manhattan Plot for Ch1 by Population")
# 
# p
# 
# 




# Done: Manhattan with ONE y-axes log scale y-axes -----------

#### My theme 
my_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(), 
    plot.subtitle = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6), # element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(.5, "mm"),
    legend.position = 'none',
    panel.border = element_blank()
  )

#### Plot function
plot_manhattan <- function(data, cutoff_data, title = "") {
  ggplot(data) +
    geom_point(
      aes(x = Cum_pos, y = log10(Score), color = Color_group, shape = Method),
      alpha = .7, size = .2
    ) +
    scale_shape_manual(
      values = c("OmegaPlus" = 16, "RAiSD" = 17)  ## Shape for Mathod
    ) +
    geom_hline(
      data = cutoff_data,
      aes(yintercept = log10(raisd), color = Population),
      linetype = "dashed", alpha = 1, linewidth = 0.2  ## Dashed cutoff for RAiSD
    ) +
    geom_hline(
      data = cutoff_data,
      aes(yintercept = log10(omegaplus), color = Population),
      linetype = "solid", alpha = 1, linewidth = 0.2  ## Solid cutoff for OmegaPlus
    ) +
    scale_x_continuous(
      # name = "Chromosome", 
      breaks = chr_label$midpoints,
      labels = chr_label$Chr,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      name = "log10(Mu or Omega)",
      expand = c(0.01, 0.01) # ,
      #sec.axis = dup_axis(
      #  # transform = ~ . / max(data$Score) * 12,  # Example transformation function
      #  name = "Omega")
    ) +
    scale_color_manual(
      values = my_colors_with_grey,
      breaks = names(my_colors),
      guide = guide_legend(override.aes = list(shape = 15, size = 1))
    ) +
    geom_vline(
      xintercept = unique(c(chr_label$start, chr_label$end)), 
      color = "grey50", linetype = "solid", alpha = 0.3
    ) +
    my_theme #+
    #ggtitle(title)
}

#### Plot p1, p2
## subset data
cultivar_list <- c("comm_Dessert", "comm_Perry")
non_cultivar_list <- c("cauc", "pyra")

if (!all(cultivar_list %in% unique(combined_data$Population))) {
  cultivar_list <- c("White", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP")
  non_cultivar_list <- c("ussu", "betu")
}

sub_combined_data <- combined_data %>% 
  filter(Population %in% cultivar_list)

sub_cutoff <- cutoff_df %>% 
  filter(Population %in% cultivar_list)

non_cultivar_data <- combined_data %>% 
  filter(Population %in% non_cultivar_list) %>%
  mutate(Color_group = ifelse(Score > Cutoff, Population, "Below_Cutoff"))

non_cultivar_cutoff <- cutoff_df  %>% 
  filter(Population %in% non_cultivar_list)

## Plot
print("Start plotting...")
p1 <- plot_manhattan(sub_combined_data, sub_cutoff, "Cultivar Populations")
p2 <- plot_manhattan(non_cultivar_data, non_cultivar_cutoff, "Non-Cultivar Populations")


combined_plot <- plot_grid(p1, p2, ncol = 1)

combined_plot
# Save
ggsave("s07.plot_Manhattan.style2/combined_plot.pdf", device = cairo_pdf, plot = combined_plot, width = 180, height = 60, units = "mm")
ggsave("s07.plot_Manhattan.style2/combined_plot.png", plot = combined_plot, width = 180, height = 60, units = "mm", dpi = 600)
ggsave("s07.plot_Manhattan.style2/combined_plot.A4.pdf", device = cairo_pdf, plot = combined_plot, width = 297, height = 210, units = "mm")
ggsave("s07.plot_Manhattan.style2/combined_plot.A4.png", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 600)
ggsave("s07.plot_Manhattan.style2/combined_plot.long.pdf", device = cairo_pdf, plot = combined_plot, width = 900, height = 60, units = "mm")
ggsave("s07.plot_Manhattan.style2/combined_plot.long.png", plot = combined_plot, width = 900, height = 60, units = "mm", dpi = 600)

 



# 
# # Delay: Manhattan with ONE y-axes only OmegaPlus, plot_grid by Pop, plot by Chr ------
# 
# #### My theme 
# my_theme <- theme_minimal() +
#   theme(
#     plot.title = element_blank(), 
#     plot.subtitle = element_blank(), 
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 6),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#     axis.text.y = element_text(size = 5),
#     legend.title = element_blank(),#element_text(size = 6),
#     legend.text = element_blank(), #element_text(size = 6),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.ticks.x = element_line(),
#     axis.ticks.length = unit(.5, "mm"),
#     legend.position = "top",
#     panel.border = element_blank()
#   )
# 
# #### Plot function
# plot_manhattan <- function(data, cutoff_data, title = "") {
#   data$Score <- log10(data$Score)                               # Log
#   cutoff_data$omegaplus <- log10(cutoff_data$omegaplus)         # Log
#   
#   ggplot(data) +
#     geom_point(
#       aes(x = Pos, y = Score, color = Color_group, shape = Method),
#       alpha = .8, size = .5
#     ) +
#     #geom_hline(
#     #  data = cutoff_data,
#     #  aes(yintercept = omegaplus, color = Population),
#     #  linetype = "solid", alpha = 1, linewidth = 0.2
#     # ) +
#     scale_x_continuous(
#       breaks = chr_label$midpoints[chr_label$Chr == unique(data$Chr)],
#       labels = chr_label$Chr[chr_label$Chr == unique(data$Chr)],
#       expand = c(0.01, 0.01)
#     ) +
#     scale_y_continuous(
#       name = "log10(Omega)",
#       expand = c(0.01, 0.01)#,
#     #  sec.axis = dup_axis(
#     #    name = "log(Omega)",
#     #    trans = "log"
#     #  )
#     ) +
#     scale_color_manual(
#       values = my_colors_with_grey,
#       breaks = names(my_colors),
#       guide = guide_legend(override.aes = list(shape = 15, size = 1))
#     ) +
#     # geom_vline(
#     #   xintercept = c(chr_label$start[chr_label$Chr == unique(data$Chr)][1], chr_label$end[chr_label$Chr == unique(data$Chr)]), 
#     #   color = "grey50", linetype = "solid", alpha = 0.5
#     # ) +
#     my_theme +
#     ggtitle(title)
# }
# 
# # For
# for (chr in paste0("Chr", 1:17)) {
#   # Select Omegaplus
#   chr_data <- combined_data %>%
#     filter(Chr == chr, Method == "OmegaPlus")
#
#
#   chr_cutoff <- cutoff_df %>%
#     filter(Population %in% unique(chr_data$Population)) %>%
#     select(Population, omegaplus)
#
#   # Check if blank
#   if (nrow(chr_data) > 0) {
#     # Split by poulation
#     plots <- list()
#     for (pop in unique(chr_data$Population)) {
#       pop_data <- chr_data %>% filter(Population == pop)
#       pop_cutoff <- chr_cutoff %>% filter(Population == pop)
#       plot_title <- paste("Manhattan Plot for", chr, " - ", pop)
#       plots[[pop]] <- plot_manhattan(pop_data, pop_cutoff, plot_title)
#     }
#     
#     # Use plot_grid arrange plot
#     combined_plot <- plot_grid(plotlist = plots, ncol = 1)
#     
#     # Save plot
#     filename <- paste0(chr, "_manhattan_plot.png")
#     ggsave(filename, plot = combined_plot, width = 16, height = 2* length(plots), dpi = 300)
#     
#     cat("Saved plot for", chr, "as", filename, "\n")
#   } else {
#     cat("No data available for", chr, ". Skipping...\n")
#   }
# }
# 
