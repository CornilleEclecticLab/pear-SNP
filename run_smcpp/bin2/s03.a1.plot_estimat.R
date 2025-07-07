# Load necessary libraries
library(ggplot2)
library(dplyr)


# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_path)

rm(list=ls())


# Read SMC++ result for plot
# data <- read.table("../output/s03.plot_estimate_cubic/s03.plot_estimate.joint.csv", header = TRUE, sep = ",")
data <- read.table("../output/s03.plot_estimate_piecewise/s03.plot_estimate.joint.csv", header = TRUE, sep = ",")
data <- data %>%
  filter(x > 0, y > 0)

# Replace 0 with 0.01 in the x column
data$x[data$x == 0] <- 0.01
data$x[data$y == 0] <- 0.01


# Read color config file
color_config <- read.table("../input/population_color.tsv", header = T)
color_uniq <- color_config %>% distinct(Population, Color)
my_color <- setNames(color_uniq$Color, color_uniq$Population)
my_color


Eu_data <- data %>% filter(label %in% c("comm_Dessert", "comm_Perry", "cauc", "pyra"))
As_data <- data %>% filter(label %in% c("betu","pash","pyri_JP","Sand_CN","Sand_CN-SE","Sand_CN-SW","ussu","White"))

# Plot curves
ggplot(As_data, aes(x = x, y = y, color = factor(label))) +
    #geom_point() +
    # geom_step() +
    geom_line() +
    # work
    # scale_x_log10() +
    # scale_y_log10() +

   # log-10 transformation introduced infinite values.
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_cartesian(xlim = c(1e3, 1e6), ylim=c(1e2,1e6)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks()+
  
  # Xilong
  # scale_x_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)), limits = c(1e3, 1e7)) +
  # scale_y_continuous(trans='log10',labels = trans_format("log10", math_format(10^.x))) +

  # Incorrect ruler
  # scale_x_continuous(trans = 'log10', labels = scales::label_math(10^.x), limits = c(1e3, 1e7)) +
  # scale_y_continuous(trans = 'log10', labels = scales::label_math(10^.x)) +

  annotate("rect", xmin = c(15000,135000), xmax = c(110000,194000), 
           ymin = 0, ymax = Inf, alpha = 0.2, fill = c("#999999", "#999999")
  ) +
  geom_vline(xintercept=19000, linetype="dashed", color = "#777777"
  ) +
  geom_vline(xintercept=33000, linetype="dashed", color = "#777777"
  ) +
  geom_vline(xintercept=2600000, color = "#777777"
  ) +

  
  scale_color_manual(values=my_color)+
  annotate(geom = "text" , x = 25000, y = 90000, 
           label = "Last glacial maximum", angle = 90, size = 3
  ) +
  annotate(geom = "text" , x = 2.38e6, y = 30000, 
           label = "PLEISTOCENE", angle = 90, size = 3
  ) +
  annotate(geom = "text" , x = 1.6e5, y = 1500, 
           label = "Penultimate glacial period", angle = 90, size = 3
  ) +
  annotate(geom = "text" , x = 45000, y = 2500, 
           label = "Last glacial period ", angle = 0, size = 3
  ) +
  
    labs(#title = "Historical Population size", 
         x = "Years Ago (gen=6 mu=3.9e-8)", y = "Ne", color = "label") +
    theme_minimal()








####################################
#>>>>>>>>>>> MSMC2 <<<<<<<<<<<<<<<<<
#>##################################

# Set mutation rate and generation time
mu <- 3.9e-08
gen <- 7.5

# Set the input folder and file names
# input_folder <- "../output/s04.run_MSMC2_within_population/"

input_folder <- "../../run_MSMC2_with_phased_data/output.backup_2024-05-24_p25-i20/s04.run_MSMC2_within_population/"
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
