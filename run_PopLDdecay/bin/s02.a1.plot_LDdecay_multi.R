library(dplyr)
library(readr)
library(ggplot2)


# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path
setwd(script_path)


# Clear up
rm(list = ls())


# Read LD decay files
ind_and_pop = read.csv("../input/s02.sample_tab_population.txt", header = FALSE, sep = '\t')
pops = unique(ind_and_pop$V2)


# Load Population order
pop_order = read.table("../input/s03.population_order.tree.txt", header=FALSE)$V1


## Load color config file
color_config <- read.table("../input/population_color.tsv", header = T)
color_uniq <- color_config %>% distinct(Population, Color)

ori_order_color <- setNames(color_uniq$Color, color_uniq$Pop)
my_color <- ori_order_color[pop_order]


# Loop through each population to load and process LD decay data
ld_data_list <- list()
for (pop in pops) {
  ld_file_path <- paste0("../output/s02.plot_LDdecay_multi/s02.plot_LDdecay_multi_r2.", pop)
  
  if (file.exists(ld_file_path)) {
    # Read LD data file
    df <- read.table(ld_file_path, header = FALSE)
    
    # Ensure correct data structure
    if (ncol(df) >= 3) {
      # Add Distance (in Kb) and LD columns, and Population as an identifier
      ld_df <- data.frame(
        Distance = df[, 1] / 1000,  # Convert to Kb
        LD = df[, 2],    # Mean r^2, or df[,3] mean D
        Population = pop
      )
      ld_data_list[[pop]] <- ld_df
    } else {
      warning(paste("The file", ld_file_path, "does not have the required columns."))
    }
  } else {
    warning(paste("File not found:", ld_file_path))
  }
}

# Combine all population data frames into one
ld_data <- bind_rows(ld_data_list)

# Create LD decay plot using ggplot
ggplot(ld_data, aes(x = Distance, y = LD, color = Population)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = my_color) +
  labs(
    title = "LD Decay across Populations",
    x = "Distance (Kb)",
    y = "LD (r^2)"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  xlim(0, 200) +
  ylim(0, 0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 200)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6))

# Save plot to PDF
pdf("../output/s02.plot_LDdecay_multi/s02.plot_LDdecay_multi_r2.colored.pdf")
print(last_plot())
dev.off()