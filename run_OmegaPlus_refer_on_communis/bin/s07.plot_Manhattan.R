library(ggplot2)
library(dplyr)
library(cowplot)


# # Set working directory to the location of the script
# script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
# script_path
# setwd(script_path)


# Clean up
rm(list=ls())

# Function to load result
load_method_data <- function(filepath, method) {
    dt <- read.table(filepath, header = FALSE)
    names(dt) <- c("Chr", "Pos", "Score")
    dt$Chr <- factor(dt$Chr, levels = paste0("Chr", 1:17))
    dt$Method <- method
    return(dt)
}


# # Load result
# omegaplus_data <- load_method_data("../input/s07.plot_Manhattan/OmegaPlus.test10000line.comm_Dessert.txt", "OmegaPlus")
# raisd_data <- load_method_data("../input/s07.plot_Manhattan/RAiSD.test10000line.comm_Dessert.txt", "RAiSD")
# sweed_data <- load_method_data("../input/s07.plot_Manhattan/SweeD.test10000line.comm_Dessert.txt", "SweeD")

# Load result
omegaplus_data <- load_method_data("../input/s07.plot_Manhattan/OmegaPlus.top20per.comm_Dessert.txt", "OmegaPlus")
raisd_data <- load_method_data("../input/s07.plot_Manhattan/RAiSD.top20per.comm_Dessert.txt", "RAiSD")
sweed_data <- load_method_data("../input/s07.plot_Manhattan/SweeD.top20per.comm_Dessert.txt", "SweeD")

# Combine data
all_data <- rbind(omegaplus_data, raisd_data, sweed_data)
raisd_data <- NULL
omegaplus_data <- NULL
sweed_data <- NULL

# Cumulate position
combined_data <- all_data %>%
    
  # Get chromosome size
  group_by(Chr) %>%
  summarise(chr_len = max(Pos)) %>%
  
  # Calulate cumulative position
  mutate(Tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this information to the data
  left_join(all_data, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos) %>%
  mutate(Cum_pos = Pos+Tot)

all_data <- NULL

axisdf <- combined_data %>%
  group_by(Chr) %>%
  summarise(mid = mean(Cum_pos))




# Load neutral envelop (cutoff values) for each method.
cutoff_raisd <- read.table('../input/s07.plot_Manhattan/comm_Dessert.RAiSD.FPR005.cutoff.txt', header = FALSE, skip=1 )
cutoff_raisd <- cutoff_raisd[1,]

cutoff_omegaplus <- read.table('../input/s07.plot_Manhattan/comm_Dessert.OmegaPlus.FPR005.cutoff.txt', header = FALSE, skip=1)
cutoff_omegaplus <- cutoff_omegaplus[1,]

cutoff_sweed <- read.table('../input/s07.plot_Manhattan/comm_Dessert.SweeD.FPR005.cutoff.txt', header = FALSE, skip=1)
cutoff_sweed <- cutoff_sweed[1,]


# RAiSD 
p1 <- ggplot(combined_data %>% filter(Method == "RAiSD"), 
             aes(x = Cum_pos, y = Score)) +
  geom_point(aes(color=Chr), alpha=0.8, size=0.7, shape=1) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 17 )) +
  scale_x_continuous(name = NULL) +   # Not show x-axis label to keep consistent
  scale_y_continuous(expand = c(0, 0), name = "RAiSD Score", limits = c(0, max(combined_data$Score[combined_data$Method == "RAiSD"]))) +
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  geom_hline(yintercept = cutoff_raisd, color = "black", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # Hide X-axis
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  )

# OmegaPlus
p2 <- ggplot(combined_data %>% filter(Method == "OmegaPlus"), 
             aes(x = Cum_pos, y = Score)) +
  geom_point(aes(color=Chr), alpha=0.8, size=0.7, shape=2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 17 )) +
  scale_x_continuous(name = NULL) +   # Not show x-axis label to keep consistent
  scale_y_continuous(expand = c(0, 0), name = "OmegaPlus Score", limits = c(0, max(combined_data$Score[combined_data$Method == "OmegaPlus"]))) +

  geom_hline(yintercept = cutoff_omegaplus, color = "black", linetype = "dashed", linewidth = 1) +
  
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # Hide X-axis
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  )


# SweeD
p3 <- ggplot(combined_data %>% filter(Method == "SweeD"), 
             aes(x = Cum_pos, y = Score)) +
  geom_point(aes(color=Chr), alpha=0.8, size=0.7, shape=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 17 )) +
  scale_x_continuous( name="Position" ,label = axisdf$Chr, breaks= axisdf$mid ) +
  scale_y_continuous(expand = c(0, 0), name = "SweeD Score", limits = c(0, max(combined_data$Score[combined_data$Method == "SweeD"]))) +

  geom_hline(yintercept = cutoff_sweed, color = "black", linetype = "dashed", linewidth = 1) +
  
  theme_minimal()+
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  )

# Combine the plots
combined_plot <- plot_grid(p1, p2, p3, ncol = 1, align = "v")
combined_plot
ggsave("combined_plot.pdf", plot = combined_plot, 
       width = 297, height = 210, units = "mm",   # horizon A4 
       device = "pdf")

q()

# OK following ###############
ggplot(combined_data, aes(x=Cum_pos, y=log10(Score)))+
  # Show all points
  geom_point(aes(color=Chr), alpha=0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # x axis:
  scale_x_continuous(label=axisdf$Chr, breaks=axisdf$mid) +
  scale_y_continuous(expand = c(0,0)) +  # Revmove space between axis and plot
  
  geom_hline(yintercept = log10(cutoff_omegaplus), color = "red", linetype = "dashed", linewidth = 1) +

  # Theme
  theme_bw() +
  theme(
  legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  
  # Add cutoff lines for each method
  #geom_hline(yintercept = log10(cutoff_raisd), color = "blue", linetype = "dashed", size = 1) +
  
  #geom_hline(yintercept = log10(cutoff_sweed), color = "red", linetype = "dashed", size = 1) +
    
  )



# OK following ###############
ggplot(combined_data, aes(x=Cum_pos, y=log10(Score))) +

  # Show all points
  geom_point( aes(color=Chr), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$mid ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

  # Custom the theme:
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

##################




# Define the shapes
shapes <- c("RAiSD" = 21, "OmegaPlus" = 22, "SweeD" = 23)

# Plot
ggplot(all_data, aes(x = cumulative_Pos, y = Score, shape = Method)) +
    geom_point(size = 2) +
    scale_shape_manual(values = shapes) +
    scale_x_continuous(breaks = chr_lengths$Chr) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "Chromosome", y = "Score", title = "Manhattan plot")


