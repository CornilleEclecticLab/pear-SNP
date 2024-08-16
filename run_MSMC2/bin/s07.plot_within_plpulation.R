library(ggplot2)
library(tidyr)
library(dplyr)

rm(list=ls())

# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_path)

# Set mutation rate and generation time
mu <- 3.9e-08
gen <- 7.5

# Set the input folder and file names
# input_folder <- "../output/s04.run_MSMC2_within_population/"
####################################
#>>>>>>>>>>> TEST <<<<<<<<<<<<<<<<<<
#>##################################
input_folder <- "../output.backup_2024-05-24_p25-i20/s04.run_MSMC2_within_population/"
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
ggplot(combined_data, aes(x = time, y = ne, color = pop, linetype = factor(rep))) +
  scale_color_manual(values=my_color)+
  geom_step() +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(x = "Years ago", y = "Effective population size", color = "Population", linetype = "Replicate") +
  theme_minimal() +
  ggtitle("MSMC2 Results: -p 1*2+25*1+1*2+1*3") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot if needed
# ggsave("msmc2_results_plot.png", width = 10, height = 7)



# Filter the data to only include population "pash"
rep1_data <- subset(combined_data, rep==1)
ggplot(rep1_data, aes(x=time, y=ne, color=pop,linetype=factor(rep)))+
  scale_color_manual(values=my_color)+
  geom_step()+
  scale_x_log10()+
  scale_y_continuous(limits=c(0, 1000)) +
  labs(x = "Years ago", y = "Effective population size", color = "Population", linetype = "Replicate") +
  theme_minimal() +
  ggtitle("MSMC2 Results: Effective Population Size Over Time") +
  theme(plot.title = element_text(hjust = 0.5))

  

# Plot Western
Eu_data <- combined_data %>% filter(pop %in% c("comm_Dessert", "comm_Perry", "cauc", "pyra"))
rep_1Eu_data <- subset(Eu_data, rep==1)
ggplot(Eu_data, aes(x = time, y = ne, color = pop, linetype = factor(rep))) +
  scale_color_manual(values=my_color)+
  geom_step() +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(x = "Years ago", y = "Effective population size", color = "Population", linetype = "Replicate") +
  theme_minimal() +
  ggtitle("MSMC2 Results: Effective Population Size Over Time") +
  theme(plot.title = element_text(hjust = 0.5))


As_data <- combined_data %>% filter(pop %in% c("betu","pash","pyri_JP","Sand_CN","Sand_CN-SE","Sand_CN-SW","ussu","White"))
rep_1As_data <- subset(As_data, rep==1)
ggplot(As_data, aes(x = time, y = ne, color = pop, linetype = factor(rep))) +
  scale_color_manual(values=my_color)+
  geom_step() +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(x = "Years ago", y = "Effective population size", color = "Population", linetype = "Replicate") +
  theme_minimal() +
  ggtitle("MSMC2 Results: Effective Population Size Over Time") +
  theme(plot.title = element_text(hjust = 0.5))







afrDat<-read.table("results/AFR.msmc2.final.txt", header=TRUE)
eurDat<-read.table("results/EUR.msmc2.final.txt", header=TRUE)
plot(afrDat$left_time_boundary/mu*gen, (1/afrDat$lambda)/(2*mu), log="x",ylim=c(0,100000),
     type="n", xlab="Years ago", ylab="effective population size")
lines(afrDat$left_time_boundary/mu*gen, (1/afrDat$lambda)/(2*mu), type="s", col="red")
lines(eurDat$left_time_boundary/mu*gen, (1/eurDat$lambda)/(2*mu), type="s", col="blue")
legend("topright",legend=c("African", "European"), col=c("red", "blue"), lty=c(1,1))