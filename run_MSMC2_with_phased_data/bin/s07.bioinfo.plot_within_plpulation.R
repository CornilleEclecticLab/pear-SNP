library(ggplot2)
library(tidyr)
library(dplyr)

rm(list=ls())

# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path
setwd(script_path)


# Set mutation rate and generation time
mu <- 3.9e-8
mu
gen <- 7.5

# Set the input folder and file names
input_folder <- "../output/s04.bioinfo.run_MSMC2_within_population/"
####################################
#>>>>>>>>>>> TEST <<<<<<<<<<<<<<<<<<
#>##################################
#input_folder <- "../output.backup_2024-05-24_p25-i20/s04.run_MSMC2_within_population/"
input_file_names <- list.files(path = input_folder, pattern = "*.final.txt", full.names = TRUE)

input_file_names

color_config <- read.table("../input/population_color.tsv", header = T)
color_uniq <- color_config %>% distinct(Population, Color)
my_color <- setNames(color_uniq$Color, color_uniq$Population)
my_color

# Initialize a list to store the data
data_list <- list()

# Read population name and replicate number from file names and populate the data list
for (i in seq_along(input_file_names)) {
  filename <- basename(input_file_names[i])
  df <- read.table(input_file_names[i], header = TRUE)
  
  # Extract population and replicate information from the filename
  split_name <- strsplit(filename, "\\.")[[1]]
  pop <- split_name[1]
  rep <- gsub("REP", "", split_name[2])
  
  # Add population and replicate information to the dataframe
  df$pop <- pop
  df$rep <- as.numeric(rep)
  
  # Transform the time and effective population size
  df$time <- df$left_time_boundary / mu * gen
  df$ne <- (1 / df$lambda) / (2 * mu)
  
  # Store the dataframe in the list
  data_list[[i]] <- df
}

# Combine all dataframes into one
combined_data <- do.call(rbind, data_list)

# Plot
# ggplot(combined_data, aes(x = time, y = ne, color = pop, linetype = factor(rep))) +
#   scale_color_manual(values=my_color)+
#   geom_step() +
#   # scale_x_continuous(limits=c(1e3,2e4))+
#   scale_x_log10() +
#   scale_y_continuous(limits = c(0, 1e5)) +
#   labs(x = "Years ago", y = "Effective population size", color = "Population", linetype = "Replicate") +
#   theme_minimal() +
#   ggtitle("MSMC2 Results:") +
#   theme(plot.title = element_text(hjust = 0.5))

# Save the plot if needed
# ggsave("msmc2_results_plot.png", width = 10, height = 7)


rep1_data <- subset(combined_data, rep==1)

# Plot Western and Eastern
# Eu_data <- combined_data %>% filter(pop %in% c("comm_Dessert", "comm_Perry", "cauc", "pyra"))
# rep1Eu_data <- subset(Eu_data, rep==1)
# As_data <- combined_data %>% filter(pop %in% c("betu","pash","pyri_JP","Sand_CN","Sand_CN-SE","Sand_CN-SW","ussu","White"))
# rep1As_data <- subset(As_data, rep==1)

ggplot(rep1_data, aes(x = time, y = ne, color = pop, linetype = factor(rep))) +
  annotate(geom = "text" , x = 25000, y = 30000, 
           label = "Last glacial maximum", angle = 90, size = 3
           ) +
  annotate(geom = "text" , x = 2.38e6, y = 30000, 
           label = "PLEISTOCENE", angle = 90, size = 3
           ) +
  annotate(geom = "text" , x = 1.4e5, y = 50000, 
           label = "Penultimate glacial period", angle = 0, size = 3
           ) +
  annotate(geom = "text" , x = 80000, y = 40000, 
           label = "Last glacial period ", angle = 0, size = 3
           ) +
  annotate("rect", xmin = c(15000,135000), xmax = c(110000,194000), 
           ymin = 0, ymax = Inf, alpha = 0.2, fill = c("#999999", "#999999")
           ) +
  geom_vline(xintercept=19000, linetype="dashed", color = "#777777"
             ) +
  geom_vline(xintercept=33000, linetype="dashed", color = "#777777"
             ) +
  geom_vline(xintercept=2600000, color = "#777777"
             ) +
  annotation_logticks()+
  scale_color_manual(values=my_color)+
  geom_step() +
  scale_x_log10() +
  scale_y_log10()+
  coord_cartesian(xlim = c(100, 3e6))+
  # scale_x_continuous(limits=c(100,1e7))+
  # scale_y_continuous(limits = c(0, 1e5)) +
  labs(x = paste0("Years ago (gen=",gen," mu=",mu,")"), y = "Effective population size", color = "Population", linetype = "Replicate") +
  theme_minimal() +
  ggtitle("MSMC2 Results: Full phased SNPs dataset") +
  theme(plot.title = element_text(hjust = 0.5))

