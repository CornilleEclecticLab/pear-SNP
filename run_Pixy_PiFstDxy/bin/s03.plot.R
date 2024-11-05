# Reference: https://pixy.readthedocs.io/en/latest/plotting.html
# Reference: https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/


## INSTALL PACKAGES ----------------------------------------------
## install CRAN packages if needed
pckgs <- c("ggplot2", "dplyr", "systemfonts", "ggforce", 
           "ggdist", "ggbeeswarm", "devtools", "ggtext")
new_pckgs <- pckgs[!(pckgs %in% installed.packages()[,"Package"])]
if(length(new_pckgs)) install.packages(new_pckgs)


## LOAD PACKAGES -------------------------------------------------
library(ggplot2) ## plotting
library(systemfonts) ## custom fonts
library(dplyr) ## data handling
# library(tidyverse)
library(readr)
library(tidyr)
library(ggdist)
library(ggtext)

## Clean up
rm(list=ls())


## CUSTOM THEME --------------------------------------------------
## overwrite default ggplot2 theme
theme_set(
  theme_minimal(
    ## increase size of all text elements
    base_size = 10, 
    ## set custom font family for all text elements
    # base_family = "Roboto Condensed"
    )
)

## overwrite other defaults of theme_minimal()
theme_update(
  ## remove major horizontal grid lines
  panel.grid.major.x = element_blank(),
  ## remove all minor grid lines
  panel.grid.minor = element_blank(),
  ## remove axis titles element_blank(), or not element_text()
  axis.title = element_text(), 
  ## larger axis text for x
  axis.text.x = element_text(size = 10),
  ## add some white space around the plot
  plot.margin = margin(rep(8, 4))
)


# Set working directory to the location of the script
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path
setwd(script_path)


# Function to convert Pixy output to long format
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  
  for (i in 1:length(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    df <- read_delim(pixy_files[i], delim = "\t")
    
    if (stat_file_type == "pi") {
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
    } else {
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
    }
    
    pixy_df[[i]] <- df
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}


# Read Pixy files from the output folder
# pixy_folder <- "../output/test.plot"
pixy_folder <- "../output/s02.run_pixy_by_chromosome.pixy_bed_1based"


pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)


## ... Waiting to load pixy data ...


# Check the total comparisons for ussu
# ussu_comp <- pixy_df %>% filter(statistic == "count_comparisons" , pop1=="ussu", is.na(pop2))
# ussu_comp$value[ussu_comp$value == "NA"] <- NA
# ussu_comp$value <- as.numeric(ussu_comp$value)
# sum_value <- sum(ussu_comp$value, na.rm = TRUE)
# print(sum_value)


# Plot: Pi ---------------------------------------------  
## Got avg_pi data by windows
pop_order = read.table("../input/s03.population_order.tree.txt", header=FALSE)$V1
pi_data <- pixy_df %>% 
  filter(statistic == "avg_pi") %>%
  mutate(pop1 = factor(pop1, levels = pop_order))
pi_data <- pi_data %>%
  filter(!is.na(value) & value != "NA")


## Calculate the final avg Pi across windows and chromosomes
diffs_data <-  pixy_df %>% 
  filter(statistic == "count_diffs") %>%
  filter(is.na(pop2)) %>%
  mutate(pop1 = factor(pop1, levels = pop_order))

sum_diffs_data <- diffs_data %>%
  group_by(pop1) %>% 
  summarise(total_value = sum(value, na.rm=TRUE)) %>% 
  mutate(total_value = as.numeric(total_value))


comparisions_data <- pixy_df %>% 
  filter(statistic == "count_comparisons") %>%
  filter(is.na(pop2)) %>%
  mutate(pop1 = factor(pop1, levels = pop_order))

sum_comparisions_data <- comparisions_data %>%
  group_by(pop1) %>% 
  summarise(total_value = sum(value, na.rm=TRUE)) %>% 
  mutate(total_value = as.numeric(total_value))


total_pi = sum_diffs_data %>% 
  inner_join(sum_comparisions_data, by="pop1", suffix = c("_diff", "_comp")) %>% 
  mutate (pi = total_value_diff / total_value_comp) %>% 
  mutate(pi = round(pi,6)) %>% 
  select(pop1, pi)


# # Create a subset to reduce file size
# n <- nrow(pi_data)
# subset_df <- pi_data[sample(n, n/20), ]


## Load color config file
color_config <- read.table("../input/population_color.tsv", header = T)
color_uniq <- color_config %>% 
  distinct(Population, Color)
  
ori_order_color <- setNames(color_uniq$Color, color_uniq$Pop)
my_color <- ori_order_color[pop_order]


## Read number of individuals (N)
Ns = read.table(
  "../output/s03.calculate_Pi_Dxy_across_chromosome.pixy_bed_1based/Pi.txt", 
  header=T)
Ns <- Ns %>% select(pop, N)
colnames(Ns) <- c("pop1", "N")


# Raincloud plot
plot <- ggplot(pi_data, aes(x=pop1, y=value, color=factor(pop1))) +
  ## half-vilolin
  ggdist::stat_halfeye(
    aes(fill = factor(pop1)),
    ## bandwidth
    adjust = .5,
    width = .6,
    ## move geom to right
    justification = -.2,
    ## remove slab interval
    .width = 0,
    point_colour = NA 
  ) +
  
  geom_boxplot(
    width = .15,
    outlier.shape = NA
  ) +

  geom_jitter(
    size = .3,
    alpha = .1,
    # control randomness and range of jitter
    position = position_jitter(seed = 0, width =.04)
  ) +
  
  # Add text for N values
  geom_text(
    data = Ns,
    aes(x = pop1, y = max(pi_data$value) * 0.75, 
        label = paste("italic(N) == ", N)),
    position = position_dodge(width = 0.8),
    parse = TRUE,
    size = 3,
    hjust = 0.5,
    vjust = 0
  ) +
  
  # Add text for total pi
  geom_text(
    data = total_pi,
    aes(x = pop1, y = pi + 0.01, 
        label = paste("pi ==", pi)),
    position = position_dodge(width = 0.8),
    parse = TRUE,
    size = 3,
    hjust = 0.5,
    vjust = 0
  ) +
  
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = .8, color = my_color),
        axis.text.y = element_text(size = 10)
  )+
  
  labs(x = "Populations", y = "Average Pi by Windows") +
  scale_color_manual(values = my_color)+
  scale_fill_manual(values = my_color) 

plot

ggsave("../output/pi_avg_boxplot.pdf", p, width = 5, height = 4)


# Statistic:  ---------------------------------------------  
# Distribution test
# shapiro.test
library(nortest)

# sample size must be between 3 and 5000 for shapiro.test() 
by(pi_data$value, pi_data$pop1, shapiro.test)
# Results: all p-values < .05, NOT following normality distribution

# use ad test
by(pi_data$value, pi_data$pop1, ad.test)
# Results: all p-values < .05, NOT following normality distribution

# use wilcox.test
library(multcompView)

tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}

pp <- pairwise.wilcox.test(pi_data$value, pi_data$pop1, p.adjust.method ="holm" )
mymat <-tri.to.squ(pp$p.value)
mymat

myletters <- multcompLetters(mymat,compare="<=",threshold=0.01,Letters=letters)

myletters_df <- data.frame(pop=names(myletters$Letters),letter = myletters$Letters )


p_sig <- plot + geom_text(data = myletters_df, mapping=aes(label = letter, x =pop, y = 0.025 ), colour="black", size=4)
p_sig
require(cowplot)

means <- aggregate(value ~ pop1, data = pi_data, FUN = mean)
sds <- aggregate(value ~ pop1, data = pi_data, FUN = sd)
result <- data.frame(treatment = means$pop1, mean = means$value, sd = sds$value)

# Box plot

plot <- ggplot() +
  geom_jitter(pi_data, mapping = aes(x = pop1, y = value), width = 0.2, color = "#DDDDDD", size = 1) +
  geom_boxplot(pi_data, mapping = aes(x = pop1, y = value, fill = pop1), outlier.shape = NA) +
  stat_summary(pi_data, mapping = aes(x = pop1, y = value, fill = pop1), fun = mean, geom = "point", shape = 23, size = 4) +
  #scale_fill_manual(values = c("#CC79A7", "#E69F00", "#D55E00", "#0072B2", "#F0E442", "#56B4E9")) +
  theme_classic() +
  ylim(0, 0.03) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1))+
  labs(x = "Populations", y = "Average Pi") 
plot

# Save the plot














######>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Convert pixy output to long format
pixy_to_long <- function(pixy_files){
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    
    } else{
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_folder <- "../output/s02.run_pixy_by_chromosome.mask_input"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)



# Plot
# Pi box plot





# Single chromosome view
# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                             avg_dxy = "D[XY]",
                             avg_wc_fst = "F[ST]"),
                             default = label_parsed)

# plotting summary statistics along a single chromosome
pixy_df %>%
  filter(chromosome == "GWHBAOS00000076") %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic))+
  geom_line(size = 0.25)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Position on Chromosome (Mb)")+
  ylab("Statistic Value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_brewer(palette = "Set1")


# A genome-wide plot of summary statistics
# create a custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                             avg_dxy = "D[XY]",
                             avg_wc_fst = "F[ST]"),
                             default = label_parsed)

# plotting summary statistics across all chromosomes
pixy_df %>%
  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                 TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:17))) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Chromsome")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
