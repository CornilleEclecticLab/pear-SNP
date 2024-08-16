library(dplyr)
library(tidyr)
library(ggplot2)

# set the working directory
the_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
show(the_path)
setwd(the_path)

rm(list = ls())

# read in and format the data
group_info <- read.table("../input/ID_vcf_vs_ID_uni.txt", header = TRUE)
group_info <- group_info %>% 
  mutate(Group = case_when(
    substr(ID_uni, 10, 10) == "W" ~ "Western",
    substr(ID_uni, 10, 10) == "E" ~ "Eastern",
    TRUE ~ "Unknown"
  ))
group_info <- group_info %>% 
  select(Sample, Group)

depth_coverage_ref_comm <- read.table("../output/s02.stat_format.pear.ref_comm.Jul2024.txt", header = TRUE)
dc_comm <- depth_coverage_ref_comm %>% 
  select(Sample, Markdu_avg_depth, Markdu_breadth_coverage)
colnames(dc_comm) <- c("Sample", "Depth_comm", "Coverage_comm")

depth_coverage_ref_pyri <- read.table("../output/s02.stat_format.pear.ref_pyri.Dec2023.txt", header = TRUE)
dc_pyri <- depth_coverage_ref_pyri %>% 
  select("Sample","MappingAvgDetpth","Breadth_Coverage")
colnames(dc_pyri) <- c("Sample", "Depth_pyri", "Coverage_pyri")


# merge the two dataframes
df <- merge(dc_comm, dc_pyri,  by = "Sample")
df <- merge(df, group_info, by = "Sample")


# Transform the data to long format
df_long <- df %>%
  pivot_longer(cols = c(Depth_comm, Depth_pyri, Coverage_comm, Coverage_pyri),
               names_to = "Variable",
               values_to = "Value")


# Split the data into two rows
df_long$Row <- ifelse(df_long$Variable %in% c("Depth_comm", "Depth_pyri"), "Depth", "Coverage")
df_long$Col <- ifelse(df_long$Variable %in% c("Depth_comm", "Coverage_comm"), "Ref P.communis", "Ref P.pyrifolia")


# Plot
ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ Row + Col, scales = "free", ncol = 2) +
  #facet_wrap(~ Row + Variable, scales = "free", ncol = 2, labeller = labeller(Row = c(co = "co1", py = "py"))) +
  theme_classic() +
  labs(title = "Depth and Coverages of Western and Eastern pears mapped onto two reference genomes",
       x = "Group",
       y = "Value") +
  theme(legend.position = "none")
