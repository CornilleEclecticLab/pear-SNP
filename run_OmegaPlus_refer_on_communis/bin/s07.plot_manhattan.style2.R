library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(data.table)

rm(list=ls())

# Set working directory to the location of the script (if running interactively)
if (interactive()) {
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
}



# Load data ---------------------------------------------------------------
#### Load RAiSD data ####
raisd_dir <- "../../run_RAiSD_refer_on_communis/output/s06.mask_results"  ## The dev folder contains a sub-results
raisd_file_list <- list.files(path = raisd_dir, pattern = "^RAiSD_Report\\..*by_win.*", full.names = TRUE)

raisd_data <- raisd_file_list %>%
  lapply(function(file_path) {
    
    # Read the data
    print(paste0("Reading file: ", file_path))
    data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE)
    
    # Check the number of columns and set column names accordingly
    if (ncol(data) == 7) {
      colnames(data) <- c("Pos", "Start", "End", "VAR", "SFS", "LD", "Score" )
    } else {
      stop(paste("Unexpected number of columns in file:", file_path))
    }
    
    # Filter out lines that are empty or start with "//"
    data <- data[!(data$Pos == "" | grepl("^//", data$Pos)), ]
    data$Pos <- as.integer(data$Pos)
    
    data$Population <- strsplit(basename(file_path), "\\.")[[1]][2]
    data$Chr <- strsplit(basename(file_path), "\\.")[[1]][3]
    
    # Set the chromosome level
    data$Chr <- factor(data$Chr, levels = paste0("Chr", 1:17))
    
    # Set Method
    data$Method <- "RAiSD"
    
    return(data)
  }) %>%
  bind_rows()

raisd_data <- raisd_data %>%
  drop_na(Score) %>%
  select(Chr, Pos, Method, Population, Score)



#### Load Omega Data ####

# Load OmegaPlus results
omega_dir <- "../output/s06.mask_results/"
omega_file_list <- list.files(path = omega_dir, pattern = "^OmegaPlus_Report.", full.names = TRUE)

omega_data <- omega_file_list %>%
  lapply(function(file_path) {
    
    # Read the data
    print(paste0("Reading file: ", file_path))
    data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE)
    
    # Check the number of columns and set column names accordingly
    if (ncol(data) == 5) {
      colnames(data) <- c("Pos", "Score", "Start", "End", "Valid")
    } else if (ncol(data) == 2) {
      colnames(data) <- c("Pos", "Score")
    } else {
      stop(paste("Unexpected number of columns in file:", file_path))
    }
    
    # Filter out lines that are empty or start with "//"
    data <- data[!(data$Pos == "" | grepl("^//", data$Pos)), ]
    data$Pos <- as.integer(as.numeric(data$Pos))
    
    data$Population <- strsplit(basename(file_path), "\\.")[[1]][2]
    data$Chr <- strsplit(basename(file_path), "\\.")[[1]][3]
    
    # Set the chromosome level
    data$Chr <- factor(data$Chr, levels = paste0("Chr", 1:17))
    
    # Set Method
    data$Method <- "OmegaPlus"
    
    return(data)
  }) %>%
  bind_rows()

# Remove duplicated ranges
omega_data <- omega_data %>%
  mutate(score_start_end = paste0(Score, "-", Start, "-", End)) %>%
  filter(!duplicated(score_start_end))

# Purge data, remove invalid score
if ("Valid" %in% colnames(omega_data)) {
  omega_data <- omega_data %>%
    filter(Valid != 0, !is.na(Valid))
}

omega_data <- omega_data %>%
  drop_na(Score) %>%
  select(Chr, Pos, Method, Population, Score)


##########  TEST  #########################################################
## /!\ Randoly select 3000 observations
set.seed(123)
omega_data_sample <- omega_data %>%
  sample_n(5000)

raisd_data_sample <- raisd_data %>%
  sample_n(5000)
###########################################################################

#### Combine data ####
# combined_data_sample <-bind_rows(raisd_data_sample, omega_data_sample)
# combined_data <- combined_data_sample
combined_data <-bind_rows(raisd_data, omega_data)





# head(combined_data)
# dim(combined_data)


#### Cumulate position ####

combined_data <- combined_data %>%
  # Get chromosome size
  group_by(Chr) %>%
  summarise(chr_len = max(Pos)) %>%
  
  # Calulate cumulative position
  mutate(Tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  
  # Add this information to data
  left_join(combined_data, ., by = c("Chr" = "Chr")) %>%
  arrange(Chr, Pos) %>%
  mutate(Cum_pos = Pos+ Tot)

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


#### Prune data ####
prune_data <- function(df, method="RAiSD", window_size=10) { # Keep one mark for every window_size (3) marks
  df <- as.data.table(df)
  df_method<- df[df$Method == method, ]
  df_other <- df[df$Method != method, ]
  
  if (nrow(df_method) == 0) {return(df)}

  df_method[, window_id := ceiling(seq_len(.N) / window_size), by = Population]
  
  pruned_df <- df_method[, .SD[which.max(Score)], by = .(Population, window_id)]
  
  final_df <- rbind(pruned_df[, !"window_id"], df_other, fill = TRUE)
  
  return(final_df)
}

## Pruning
print(paste0("Pruning RAiSD data from  ", nrow(combined_data), "..."))
combined_data <- prune_data(combined_data, "RAiSD", 20)
print(paste0("Done, purning data, left ", nrow(combined_data), "..." ))

#### Load cutoff ####
populations <- na.omit(unique(combined_data$Population))

# Initialize
cutoff_df <- data.frame(Population = character(),  raisd = numeric(),  omgegaplus = numeric(),  stringsAsFactors = FALSE)

# Loop through populations and read cutoff values
for (pop in populations) {
  # Generate file paths for each method
  raisd_file <- paste0("../../run_MS_simulate_refer_on_communis/output/s03.find_FPR_cutoff_RAiSD/", pop, ".RAiSD.FPR0005.cutoff.txt")
  omegaplus_file <- paste0("../../run_MS_simulate_refer_on_communis/output/s03.find_FPR_cutoff_OmegaPlus/",pop, ".OmegaPlus.FPR0005.cutoff.txt")

  # Read cutoff values from each file (taking only the first value)
  raisd_cutoff <- read.table(raisd_file, header = FALSE, skip = 1)[1, 1]
  omegaplus_cutoff <- read.table(omegaplus_file, header = FALSE, skip = 1)[1, 1]

  # Add to the dataframe
  cutoff_df <- rbind(
    cutoff_df,
    data.frame(
      Population = pop,
      raisd = raisd_cutoff,
      omegaplus = omegaplus_cutoff
    )
  )
}


combined_data <- combined_data %>%
  left_join(cutoff_df, by = "Population") %>%
  mutate(
    Cutoff = case_when(
      Method == "RAiSD" ~ raisd,
      Method == "OmegaPlus" ~ omegaplus
    ),
    Color_group = ifelse(Score > Cutoff, Population, "Below_Cutoff"),
    Color_group = factor(Color_group, levels = c(unique(Population), "Below_Cutoff"))
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

my_colors_with_grey <- c(my_colors, "Below_Cutoff" = "#888888")


# Manhattan with dual y-axes ----------------------------------------------

#### My theme ####
my_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(), # element_text(size = 8, face = "bold"),
    plot.subtitle = element_blank(), # element_text(size = 8),
    axis.title = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(.5, "mm"),
    legend.position = "top",
    panel.border = element_blank()
  )


#### Plot function ####
plot_manhattan <- function(data, cutoff_data, title = "") {
  ggplot(data) +
    geom_point(
      aes(x = Cum_pos, y = Score, color = Color_group, shape = Method),
      alpha = .8, size = .5
    ) +
    geom_hline(
      data = cutoff_data,
      aes(yintercept = raisd, color = Population),
      linetype = "dashed", alpha = 1, linewidth = 0.2
    ) +
    geom_hline(
      data = cutoff_data,
      aes(yintercept = omegaplus, color = Population),
      linetype = "solid", alpha = 1, linewidth = 0.2
    ) +
    scale_x_continuous(
      name = "Chromosome", 
      breaks = chr_label$midpoints,
      labels = chr_label$Chr,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      name = "Mu",
      expand = c(0.01, 0.01),
      sec.axis = dup_axis(
        # transform = ~ . / max(data$Score) * 12,  # Example transformation function
        name = "Omega")
    ) +
    scale_color_manual(
      values = my_colors_with_grey,
      breaks = names(my_colors),
      guide = guide_legend(override.aes = list(shape = 15, size = 1))
    ) +
    geom_vline(
      xintercept = c(chr_label$start[1], chr_label$end), 
      color = "grey50", linetype = "solid", alpha = 0.5
    ) +
    my_theme +
    ggtitle(title)
}

#### Plot p1, p2 ####
## subset data
cultivar_list <- c("comm_Dessert", "comm_Perry")

sub_combined_data <- combined_data %>% 
  filter(Population %in% cultivar_list)

sub_cutoff <- cutoff_df %>% 
  filter(Population %in% cultivar_list)

non_cultivar_data <- combined_data %>% 
  filter(!Population %in% cultivar_list) %>%
  mutate(Color_group = ifelse(Score > Cutoff, Population, "Below_Cutoff"))

non_cultivar_cutoff <- cutoff_df  %>% 
  filter(!Population %in% cultivar_list)

## Plot
print("Start plotting...")
p1 <- plot_manhattan(sub_combined_data, sub_cutoff, "Cultivar Populations")
p2 <- plot_manhattan(non_cultivar_data, non_cultivar_cutoff, "Non-Cultivar Populations")


combined_plot <- plot_grid(p1, p2, ncol = 1)

combined_plot


# Save
ggsave("s07.plot_Manhattan.style2/combined_plot.pdf", plot = combined_plot, width = 180, height = 80, units = "mm", dpi = 300)
ggsave("s07.plot_Manhattan.style2/combined_plot.png", plot = combined_plot, width = 180, height = 80, units = "mm", dpi = 300)
ggsave("s07.plot_Manhattan.style2/combined_plot.A4.pdf", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 300)
ggsave("s07.plot_Manhattan.style2/combined_plot.A4.png", plot = combined_plot, width = 297, height = 210, units = "mm", dpi = 300)



