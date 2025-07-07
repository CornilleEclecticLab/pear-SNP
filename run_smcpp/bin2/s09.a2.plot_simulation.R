library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)

rm(list=ls())

# Function to read and process the data
create_boxplots <- function(data_file) {
  df <- read.table(data_file, header = TRUE)

  # Select the columns to plot
  cols_to_plot <- c("theat_hat_w", "pi", "ho", "he", "fa")

  # Convert data to long format for ggplot
  df_long <- df %>%
    select(population, all_of(cols_to_plot)) %>%
    gather(key = "Variable", value = "Value", -population) %>%
    # Create a flag for real data
    mutate(is_real = population == "real")
  
  # Create the boxplot with individual points
  ggplot(df_long, aes(x = Variable, y = Value)) +
    # Add regular boxplot
    # geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    # Add simulated data points
    geom_jitter(data = subset(df_long, !is_real),
                width = 0.2, alpha = 0.4, color = "blue", size = 1) +
    # Add real data points with different color and size
    geom_point(data = subset(df_long, is_real),
               color = "purple", size = 3, shape = 17) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      legend.position = "top"
    ) +
    labs(
      title = "Distribution of Genetic Variables in Simulated and Real Data",
      subtitle = "Purple triangles indicate real data points",
      x = "Variables",
      y = "Values"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

  # Save the plot
  # ggsave("genetic_variables_boxplot.pdf", width = 10, height = 6)
}

# Usage:
create_boxplots("../output/s09.dev_OK.assess_simulation_pi/msprime.comm_Dessert.REPs.norm.vcf.gz.stat.summary.csv")
