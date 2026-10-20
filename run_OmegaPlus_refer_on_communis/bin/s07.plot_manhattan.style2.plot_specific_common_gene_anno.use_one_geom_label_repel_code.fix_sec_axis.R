suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(data.table)
  # library(extrafont)
})

# If run on Linux, where Arial front is not installed, using this to add it. 
# But the text will be converted to vector objects, not characters, in the output PDF file.

# library(showtext)
# font_add("Arial", "./s07.font/arial.ttf")
# showtext_auto()


rm(list=ls())

# Set working directory to the location of the script (if running interactively)
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}


binary_cache_file <- "combined_data_Omega_TajimaD_binary.rds"
load(binary_cache_file)
print("Data loaded from binary file")

unique(combined_data$Method)

# Reduce data, and remove 0 for log
combined_data <- combined_data %>%
  filter(
    (Method == "OmegaPlus" & Score >= 1) |
    (Method == "Pixy_TajimaD" & Score <= -0.5)
  )

# Edit the Color_group column
combined_data <- combined_data %>%
  mutate(
    Chr_num = as.numeric(gsub("[^0-9]", "", Chr)),  
    Color_group = case_when(
      Color_group == "Below_Cutoff" & Chr_num %% 2 == 1 ~ "Below_Cutoff_odd",
      Color_group == "Below_Cutoff" & Chr_num %% 2 == 0 ~ "Below_Cutoff_even",
      TRUE ~ Color_group
    )
  )

#Convert negative TajimaD values to positive, and create logScore for plotting
combined_data <- combined_data %>%
  mutate(
    Score_trans = case_when(
      Method == "OmegaPlus" ~ Score,
      Method == "Pixy_TajimaD" ~ -Score  # Flip TajimaD so it appears above x-axis
    ),
    Cutoff_trans = case_when(
      Method == "OmegaPlus" ~ Cutoff,
      Method == "Pixy_TajimaD" ~ -Cutoff  # Flip TajimaD cutoff so it appears above x-axis
    ),
    logScore = log10(Score_trans),  # Avoid log10(0)
    logCutoff = log10(Cutoff_trans),
    log_omegaplus = log10(omegaplus)
  )

# Rescale for dual y-axis plotting only
combined_data$Rescaled_times <- combined_data$log_omegaplus / combined_data$logCutoff

combined_data$logCutoff <- combined_data$logCutoff * combined_data$Rescaled_times

combined_data$logScore <- combined_data$logScore * combined_data$Rescaled_times

# Drop some small data:      log10(score)<0.4
combined_data <- combined_data %>%
  mutate(Cutoff_trans = case_when(
    Method == "OmegaPlus" ~ Cutoff,
    Method == "Pixy_TajimaD" ~ -Cutoff
  )) %>%
  filter(!(logScore < 0.6 & Score_trans < Cutoff_trans))  # Remove data with logScore < 0.6 to save space in the plot


head(combined_data)

# Prepare transformed cutoff values
cutoff_df$log_tajimad <- log10(-cutoff_df$tajimad)
cutoff_df$log_omegaplus <- log10(cutoff_df$omegaplus)
cutoff_df$Rescaled_times <- cutoff_df$log_omegaplus / cutoff_df$log_tajimad
cutoff_df$logCutoff <- cutoff_df$log_omegaplus


############## TEST ################
# if (interactive()) {
#   set.seed(123)
#   sampled_data <- combined_data[sample(nrow(combined_data), 100000), ]
#   combined_data <- sampled_data
#   print("INFO: on smapled dataset")
# }

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

my_colors_with_grey <- c(my_colors, "Below_Cutoff" = "#AAAAAA",
                         "Below_Cutoff_odd" = "#AAAAAA",
                         "Below_Cutoff_even" = "grey85",
                         "OmegaPlus_sig" = "#000000")  # NEw: black outlier circles


# Color key for OmegaPlus in black
combined_data <- combined_data %>%
  mutate(
    Color_key = case_when(
      grepl("^Below_Cutoff", Color_group) ~ Color_group,
      Method == "OmegaPlus" ~ "OmegaPlus_sig",
      TRUE ~ Population
    )
  )


# Add gene anno -----------------------------------------------------------
## Read gene annotation
library(readr)
common_gene_annotations <- read_tsv(
  "../output/s08.summary_table/manually_interest_common_pop.positive_selection_genes_summary.comm_Dessert.comm_Perry.txt") %>% 
  select(Gene_ID, Chr_ID, Start, End, Common_vs_specific_in_population, Symbol) %>%
  filter(!(Symbol %in% c("-", "NA"))) %>%
  drop_na()
print(common_gene_annotations)


# Merge combined_data with gene_annotations based on matching Chr
common_gene_annotations <- common_gene_annotations %>%
  dplyr::rename(gene_start = Start, gene_end = End)

combined_with_common_genes <- combined_data %>%
  as.data.frame() %>%
  filter(case_when(
    Method == "OmegaPlus" ~ Score >= Cutoff,
    Method == "Pixy_TajimaD" ~ Score <= Cutoff
  )) %>%
  dplyr::left_join(common_gene_annotations, by = c("Chr" = "Chr_ID"), relationship = "many-to-many") %>%
  filter(Start <= gene_end & End >= gene_start) %>%
  select(-gene_start, -gene_end) %>%
  group_by(Gene_ID, Population) %>%
  filter(Score == max(Score)) %>%  # Keep the row with the highest Score for each Gene_ID and Population
  ungroup()  # Ungroup after filtering to return to a regular data frame


# View the resulting table
head(combined_with_common_genes)



# Load population specific gene annotation
specific_gene_annotations <- read.table("../output/s08.summary_table/manually_specific.positive_selection_genes.summary.txt",
                                   sep='\t', header=TRUE, quote = "")%>% 
  select(Gene_ID, Chr_ID, Start, End, Common_vs_specific_in_population, Detected_populations, Symbol) %>%
  filter(!(Symbol %in% c("-", "NA"))) %>% 
  na.omit() 

specific_gene_annotations <- specific_gene_annotations %>%
  dplyr::rename(gene_start = Start, gene_end = End)

# Merge combined_data with based on matching Chr
combined_with_specific_genes <- combined_data %>%
  as.data.frame() %>%
  filter(case_when(
    Method == "OmegaPlus" ~ Score >= Cutoff,
    Method == "Pixy_TajimaD" ~ Score <= Cutoff
  )) %>%
  dplyr::left_join(
    specific_gene_annotations,
    by = c("Chr" = "Chr_ID",
            "Population" = "Detected_populations"),
    relationship = "many-to-many") %>%
  filter(Start <= gene_end & End >= gene_start) %>%
  select(-gene_start, -gene_end) %>%
  group_by(Gene_ID, Population) %>%
  filter(Score == max(Score)) %>%  # Keep the row with the highest Score for each Gene_ID and Population
  ungroup()  # Ungroup after filtering to return to a regular data frame

head(combined_with_specific_genes)


# Combine the common and specific genes
combined_with_genes <- bind_rows(combined_with_common_genes, combined_with_specific_genes)
combined_with_genes <- combined_with_genes %>%
  mutate(
    label_fill = case_when(
      Common_vs_specific_in_population == "specific" ~ my_colors_with_grey[Population],
      Common_vs_specific_in_population == "common" ~ "white"
    ),
    label_color = case_when(
      Common_vs_specific_in_population == "specific" ~ "white",
      Common_vs_specific_in_population == "common" ~ my_colors_with_grey[Population]
    )
  )

head(combined_with_genes)



library(ggrepel)
#### My theme 
my_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(), 
    plot.subtitle = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5, margin = margin(t = 0)),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6), # element_blank(), 
    text = element_text(size = 6, family = "Arial"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(), #element_line(),
    # axis.ticks.length = unit(.5, "mm"),
    legend.position = 'none',
    panel.border = element_blank()
  )

#### Plot function
plot_manhattan <- function(data, cutoff_data, title = "", gene_annotaion_data) {
  # Separate significant and non-significant data
  data_significant <- data %>% filter(!grepl("^Below_Cutoff", Color_group))
  data_nonsignificant <- data %>% filter( grepl("^Below_Cutoff", Color_group))
  pop_color <- my_colors_with_grey[data_significant$Population[1]]
  rescale_factor <- cutoff_data$Rescaled_times[1]
  
  ggplot() +
    # Below cutoff
    geom_point(
      data = data_nonsignificant,
      aes(x = Cum_pos, y = logScore, color = Color_key, shape = Method),
      alpha = 0.8, size = 0.3
    ) +
    # Over cutoff 
    geom_point(
      data = data_significant,
      aes(x = Cum_pos, y = logScore, color = Color_key, shape = Method),
      alpha = .8, size = 0.5
    ) +
    # Define the shapes for scattrers
    scale_shape_manual(values = c("OmegaPlus" = 16, "Pixy_TajimaD" = 17)  ## Shape for Mathod
    ) +
    # Add cutoff lines
    geom_hline(
      data = cutoff_data,
      aes(yintercept = logCutoff, color = "grey"),
      linetype = "dashed", alpha = 1, linewidth = 0.2  ## Dashed cutoff for RAiSD
    ) +
    # geom_hline(
    #   data = cutoff_data,
    #   aes(yintercept = log10(omegaplus), color = "grey"),
    #   linetype = "solid", alpha = 1, linewidth = 0.2  ## Solid cutoff for OmegaPlus
    # ) +
    # Add cutoff line annotations
    # annotate("text", 
    #          x = max(data$Cum_pos) * 0.99, 
    #          y = log10(cutoff_data$tajimad[1]), 
    #          label = "Mu cutoff", 
    #          hjust = 1, vjust = -0.5, 
    #          size = 2, color = pop_color
    # ) +
    # annotate("text", 
    #          x = max(data$Cum_pos) * 0.99, 
    #          y = log10(cutoff_data$omegaplus[1]), 
    #          label = "Omega cutoff", 
    #          hjust = 1, vjust = -0.5, 
    #          size = 2, color = pop_color
    # ) +
    
    # Add gene annotations
    # geom_text_repel(
    geom_label_repel(
      data = gene_annotaion_data, 
      aes(x = Cum_pos, y = logScore, label = paste0("italic('", Symbol, "')")), 
      parse = TRUE,
      fill = gene_annotaion_data$label_fill,
      color = gene_annotaion_data$label_color,
      #color = pop_color, 
      size = 0.5*.pt,
      family = "Arial",
      alpha = 0.7,
      force = 0.4, #0.1,
      segment.color = pop_color,
      min.segment.length = 0.05, segment.size = 0.1, #.2,
      label.padding = 0.08, # 0.15,
      # box.padding = 0.15, # 0.1,    
      point.padding = 0.1,  
      max.iter = 3000,
      nudge_x = 0,          
      nudge_y = 0.05,       
      segment.curvature = 0,
      max.time = 2,         
      max.overlaps = Inf
    ) +
    # The chromosome boundary line
    # geom_vline(
    #   xintercept = unique(c(chr_label$start, chr_label$end)), 
    #   color = "grey", linetype = "solid", alpha = 0.3
    # ) +
    # Scales and theme
    scale_x_continuous(
      # name = "Chromosome", 
      breaks = chr_label$midpoints,
      labels = chr_label$Chr,
      expand = c(0.001, 0.001)
    ) +
    scale_y_continuous(
      # name = "log10(Omega)",
      # name = "● log₁₀(Omega)",
      # name = bquote("\u25CF"~~log[10]~"(Omega)"),     # ●
      name = bquote("\u25CF"~~log[10]~omega),
      
      
      #sec.axis = dup_axis(name = bquote("\u25B2"~~log[10]~"(-Tajima's D)")),  # ▲
      sec.axis = sec_axis(~ . / rescale_factor, name = bquote("\u25B2"~~log[10]*"("~"-Tajima's "*italic(D)*")")),
      
      # name = bquote(log[10]~"(Omega)"),
      # sec.axis = dup_axis(name = bquote(log[10]~"(Rescaled Mu)")),  # ▲

      # expand = expansion(add = c(0.01, 0.1))
      expand = expansion(mult = c(0, 0.1), add = c(0, 0))
      # expand = c(-0.1, 0.2) # ,
      # sec.axis = dup_axis(
      #  # transform = ~ . / max(data$Score) * 12,  # Example transformation function
      # name = "Omega")
    ) +
    scale_color_manual(
      values = my_colors_with_grey,
      breaks = c("OmegaPlus_sig", names(my_colors))
      # guide = guide_legend(override.aes = list(shape = 15, size = 1))
    ) +
    my_theme +
    theme(
      axis.title.y.left  = element_text(size = 6,color = "black"),
      axis.title.y.right = element_text(size = 6, color = pop_color)
    ) #+
  #ggtitle(title)
}

#### Plot p1, p2
## subset data

cultivar_p1 <- "comm_Dessert"
cultivar_p2 <- "comm_Perry"


sub_combined_data_p1 <- combined_data %>% 
  filter(Population == cultivar_p1)

sub_cutoff_p1 <- cutoff_df %>% 
  filter(Population == cultivar_p1)

sub_combined_with_gene_p1 <- combined_with_genes %>% 
  filter(Population == cultivar_p1)

sub_combined_data_p2 <- combined_data %>% 
  filter(Population == cultivar_p2)
sub_cutoff_p2 <- cutoff_df %>% 
  filter(Population == cultivar_p2)
sub_combined_with_gene_p2 <- combined_with_genes %>% 
  filter(Population == cultivar_p2)


## Plot
print("Start plotting...")
p1 <- plot_manhattan(sub_combined_data_p1, sub_cutoff_p1, cultivar_p1, sub_combined_with_gene_p1)
p1

p2 <- plot_manhattan(sub_combined_data_p2, sub_cutoff_p2, cultivar_p2, sub_combined_with_gene_p2)
p2
combined_plot <- plot_grid(p1, p2, ncol = 1)
combined_plot

#### Save plot ---------
ggsave("s07.plot_Manhattan.style2.1/combined_plot.pdf", device = cairo_pdf, plot = combined_plot, width = 130, height = 60, units = "mm")
# ggsave("s07.plot_Manhattan.style2.1/combined_plot.png", plot = combined_plot, width = 130, height = 60, units = "mm", dpi = 1200)

# ggsave("s07.plot_Manhattan.style2.1/combined_plot.A4.pdf", device = cairo_pdf, plot = combined_plot, width = 297, height = 210, units = "mm")
# ggsave("s07.plot_Manhattan.style2.1/combined_plot.A4.png", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 1200)

# ggsave("s07.plot_Manhattan.style2.1/combined_plot.long.pdf", device = cairo_pdf, plot = combined_plot, width = 900, height = 60, units = "mm")
# ggsave("s07.plot_Manhattan.style2.1/combined_plot.long.png", plot = combined_plot, width = 900, height = 60, units = "mm", dpi = 1200)

print("C'est bon!")

# Upset plot  -----------------------------------------------------------------
# Read the file list
W_file_list <- readLines("../../run_OmegaPlus_refer_on_communis/input/s08.summary_list.txt")

# Function read the data
read_data <- function(file_path) {
  print(paste0("Reading file: ", file_path))
  
  if (file.size(file_path) == 0) {
    print(paste0("File ", file_path, " is empty. Skipping..."))
    return(NULL)
  }
  
  # Read file content
  first_line <-  read.table(file_path, header = FALSE, sep = '\t', fill = TRUE, nrows = 1)
  print(paste0("First line of ", file_path, ": ", paste(first_line[1, ], collapse = "\t")))
  
  if (ncol(first_line) == 2) {
    data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, col.names = c("geneID", "Method"))
    data$Population <- strsplit(basename(file_path), "\\.")[[1]][3]
  } else {
    data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, col.names = c("geneID", "Method", "Population"))
  }
  
  return(data)
}

# Read all files and combine them
W_df <- lapply(W_file_list, read_data) %>% bind_rows()

# Load color
populations_in_data <- unique(W_df$Population)
W_my_colors <- ori_order_color[names(ori_order_color) %in% populations_in_data]
print("Loading color config")
W_my_colors



# Convert data into presence-absence matrix ------------------------------------
W_wider_df <- W_df %>%
  select(geneID, Population) %>%
  distinct() %>%
  mutate(present = 1) %>% 
  pivot_wider(names_from = Population, values_from = present, values_fill = 0)

# Convert tibble to data frame
W_wider_df <- as.data.frame(W_wider_df)

# Set geneID as row names
rownames(W_wider_df) <- W_wider_df$geneID
W_wider_df <- W_wider_df[, -1]


# Plot using complexUpset ------------------------------------------------------
library(ComplexUpset)

# Change name
W_wider_new_name_df <- W_wider_df %>%
  dplyr::rename(
    `P. communis Dessert` = `comm_Dessert`,
    `P. communis Perry`   = `comm_Perry`,
    `P. pyraster`         =  pyra,
    `P. caucasica`        =  cauc)

W_set_order <- c(
  "P. communis Dessert",
  "P. communis Perry",
  "P. pyraster",
  "P. caucasica"
)
W_set_order <- rev(W_set_order)


W_p <- ComplexUpset::upset(
  W_wider_new_name_df, 
  intersect = W_set_order,
  min_size = 4,
  set_sizes = FALSE,
  sort_set = FALSE, # set F to use custom order
  queries = list(
    upset_query(set = 'P. communis Dessert', fill = "#442CF4", color = "#442CF4"),
    upset_query(set = 'P. communis Perry',   fill = "#648FFF", color = "#648FFF"),
    upset_query(set = 'P. pyraster',         fill = "#035b03", color = "#035b03"),
    upset_query(set = 'P. caucasica',        fill = "#27bc40", color = "#27bc40")
  ),
  base_annotations = list(
    'Gene count' = intersection_size(
      text = list(
        vjust = -0.2,
        hjust = 0.1,
        angle = 45,
        size = 1.5
      )
    )
  ),
  width_ratio = 0.2,
  height_ratio = 0.7, 
  matrix = intersection_matrix(
    geom = geom_point(),#shape = 21, size = 1.5),
    outline_color = list(active = NA, inactive = 'grey'),
    segment = geom_segment(linetype='solid', color = 'black', linewidth = 0.3)
  ),
  # )+ coord_cartesian(clip = 'off'),
  themes = upset_modify_themes(
     list('intersections_matrix' = theme(axis.text.y = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.y = element_text(size=6, angle=90),
                                         panel.grid.major.x = element_blank(),
                                         panel.border = element_blank(),
                                         text = element_text(family='Arial'),
                                         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
                                         ),
          'Gene count' =  theme(axis.title.y = element_text(size = 6),
                                axis.text.y = element_text(size = 5),
                                panel.background = element_blank(),
                                panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.2),
                                panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.2),
                                panel.grid.minor.y = element_line(color = "#EBEBEB", linewidth = 0.1),
                                panel.border = element_blank(),
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks = element_blank(),
                                # text = element_text(family='Arial'),
                                plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
                              ) 
   )
  )
)+ 
  labs(y = "Populations") &
  coord_cartesian(clip = 'off') &
  theme(text = element_text(family='Arial'))
W_p


print("Start Upset plot")

combin_upset_manhat <-plot_grid(W_p, combined_plot, ncol = 2,  rel_widths = c(1,2.3))

combin_upset_manhat



#### Save plot ---------
ggsave("s07.plot_Manhattan.style2.1/combined_upset_manhattan_plot.pdf", device = cairo_pdf, plot = combin_upset_manhat, width = 178, height = 65, units = "mm")
ggsave("s07.plot_Manhattan.style2.1/combined_upset_manhattan_plot.svg", plot = combin_upset_manhat, width = 178, height = 65, units = "mm")
ggsave("s07.plot_Manhattan.style2.1/combined_upset_manhattan_plot.png", plot = combin_upset_manhat, width = 178, height = 65, units = "mm", dpi = 1200)

print("C'est bon!")
