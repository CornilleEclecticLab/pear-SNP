# By Xilong CHEN
# Create date: 2023-09-21
# Contact: chen_xilong@outlook.com

# Modified by Yuqi on 2025-06-11, 06-16

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(cowplot)

# Set working directory to the location of the script and clean up
if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}


###################
## Homo vs He
###################

### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)


#df_sweep <- df_sweep %>%  
#    filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "cauc", "pyra",
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c(
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP",
  "White", "pash", "ussu", "betu"
)

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection")))


### use all_masked regions of the dm
df_sweep_mask <- df_sweep[df_sweep$region_type == "all_masked", ]



# get max value for p position
ctrl_max <- df_sweep_mask %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05


p <- ggboxplot(df_sweep_mask, x = "pop", y = "de_nu",
               color = "type", palette = "aaas") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add p-values
p + stat_compare_means(
  aes(group = type),
  label       = "p.format",
  label.y  = ctrl_pos
) 



#### plot seperately for West and East ####
# rm All
df_sweep_mask <- df_sweep_mask %>%
  filter(type != "All")


# get max value for p position
ctrl_max_w <- df_sweep_mask %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# Plot west
p_w <- ggboxplot(
    data = df_sweep_mask %>% filter(pop %in% west_pops),
    x = "pop", 
    y = "de_nu",
    color = "type", 
    palette = "jco",
    add = "jitter"
  ) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
    annotate("text", x = 1.5, y=60, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 3.5, y=60, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Occidental populations", y = "No. of Deleterious alleles", color="Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      aes(group = type),
      label       = "p.format",
      label.y  = ctrl_pos_w)
p_w

## East
# get max value for p position
ctrl_max_e <- df_sweep_mask %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_mask %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_nu",
    color = "type", 
    palette = "jco",
    add = "jitter"
  ) +
    geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
    annotate("text", x = 2.5, y=75, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 7, y=75, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Oriental populations", y = NULL, color="Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = type),
    #comparisons = list(c("selection", "control")),
    label = "p.format",
    label.y = ctrl_pos_e
  )

p_m <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m



###################
## per CDS
###################

### read data
rm(list=ls())
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)


df_sweep <- df_sweep %>%  
  filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "cauc", "pyra",
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c(
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP",  "White", 
  "pash", "ussu", "betu"
)

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection")))


df_sweep$de_per_100k_cds <- df_sweep$de_per_cds * 100000


### use all alleles of the dm
df_sweep_all <- df_sweep[df_sweep$type == "All", ]


df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]

my_comparisons <- list( c("all_masked", "control"), c("control", "selection"), c("all_masked", "selection") )

# get max value for p position
ctrl_max <- df_sweep_all_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05


p <- ggboxplot(df_sweep_all_sel_con, x = "pop", y = "de_per_100k_cds",
               color = "region_type", palette = "npg") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles/CDS length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add p-values
p + stat_compare_means(
  aes(group = region_type),
  label       = "p.format",
  label.y  = ctrl_pos
) 



#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_all_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# Plot west
p_w <- ggboxplot(
    data = df_sweep_all_sel_con %>% filter(pop %in% west_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "region_type", 
    palette = "npg",
    add = "jitter"
  ) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
    annotate("text", x = 1.5, y=60, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 3.5, y=60, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/CDS length", color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      aes(group = region_type),
      label       = "p.format",
      label.y  = ctrl_pos_w)
p_w

## East
# get max value for p position
ctrl_max_e <- df_sweep_all_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_all_sel_con %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "region_type", 
    palette = "npg",
    add = "jitter"
  ) +
    geom_vline(xintercept = 4.5, linetype = "dashed", color = "black") +
    annotate("text", x = 2.5, y=75, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 5.5, y=75, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Oriental populations", y = NULL, color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.format",
    label.y = ctrl_pos_e
  )

p_m <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m

ggsave2("../output/s06.dm_sweep/pop_allele_CDS.region.mutation_burden.svg", p_m, width=10, height=3)
ggsave2("../output/s06.dm_sweep/pop_allele_CDS.region.mutation_burden.pdf", p_m, width=10, height=3)









### plot only Dessert and Perry
df_sweep_all_dom <- df_sweep_all[df_sweep_all$pop == "comm_Dessert" | df_sweep_all$pop == "comm_Perry", ]

df_sweep_all_dom$pop <- factor(df_sweep_all_dom$pop, c("comm_Dessert", "comm_Perry"))

df_sweep_all_dom$region_type <- factor(df_sweep_all_dom$region_type, c("all_masked", "control", "selection"))

df_sweep_all_dom$de_per_100k_cds <- df_sweep_all_dom$de_per_cds * 100000

p <- ggboxplot(df_sweep_all_dom, x = "pop", y = "de_per_100k_cds",
               color = "region_type", palette = "npg"
) + theme_bw(base_size = 8)
# Use only p.format as label. Remove method name.
p + stat_compare_means(aes(group = region_type), label="p.format" )#, label = "p.signif")

ggsave2("../output/s06.dm_sweep/pop_comm.region.mutation_burden.svg", width=5, height=2.5)
ggsave2("../output/s06.dm_sweep/pop_comm.region.mutation_burden.pdf", width=5, height=2.5)










###############
## per SNP ----
###############

rm(list=ls())

### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.SNP.pop_mutation_burden.Rdata.txt", header =T)


df_sweep <- df_sweep %>%  
    filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "cauc", "pyra",
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c(
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP",
  "White", "pash", "ussu", "betu"
)

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection")))


### use all alleles of the dm
df_sweep_all <- df_sweep[df_sweep$type == "All", ]

df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]



# get max value for p position
ctrl_max <- df_sweep_all_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_snp, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05


p <- ggboxplot(df_sweep_all_sel_con, x = "pop", y = "de_per_snp",
               color = "region_type", palette = "npg") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles/SNV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add p-values
p + stat_compare_means(
  aes(group = region_type),
  label       = "p.format",
  label.y  = ctrl_pos
) 



#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_all_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_snp, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# Plot west
p_w <- ggboxplot(
    data = df_sweep_all_sel_con %>% filter(pop %in% west_pops),
    x = "pop", 
    y = "de_per_snp",
    color = "region_type", 
    palette = "npg",
    add = "jitter"
  ) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
    annotate("text", x = 1.5, y=0.25, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 3.5, y=0.25, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/SNVs", color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(
      aes(group = region_type),
      label       = "p.format",
      label.y  = ctrl_pos_w)
p_w

## East
# get max value for p position
ctrl_max_e <- df_sweep_all_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_snp, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_all_sel_con %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_snp",
    color = "region_type", 
    palette = "npg",
    add = "jitter"
  ) +
    geom_vline(xintercept = 4.5, linetype = "dashed", color = "black") +
    annotate("text", x = 2.5, y=0.28, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 5.5, y=0.28, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Oriental populations", y = NULL, color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.format",
    label.y = ctrl_pos_e
  )
p_e

p_m <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m

ggsave2("../output/s06.dm_sweep/pop_allele_SNP.region.mutation_burden.svg", p_m, width=10, height=3)
ggsave2("../output/s06.dm_sweep/pop_allele_SNP.region.mutation_burden.pdf", p_m, width=10, height=3)





#################
# synonymyous SNP
#################
rm(list=ls())

### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt", header =T)


df_sweep <- df_sweep %>%  
  filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "cauc", "pyra",
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c(
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP",
  "White", "pash", "ussu", "betu"
)

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection")))


### use all alleles of the dm
df_sweep_all <- df_sweep[df_sweep$type == "All", ]

df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]


# get max value for p position
ctrl_max <- df_sweep_all_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_syno, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05


p <- ggboxplot(df_sweep_all_sel_con, x = "pop", y = "de_per_syno",
               color = "region_type", palette = "npg") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles/SNV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add p-values
p + stat_compare_means(
  aes(group = region_type),
  label       = "p.format",
  label.y  = ctrl_pos
) 


#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_all_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_syno, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# Plot west
p_w <- ggboxplot(
  data = df_sweep_all_sel_con %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_syno",
  color = "region_type", 
  palette = "npg",
  add = "jitter"
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y=0.7, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y=0.7, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/synonymyous SNVs", color="Region type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    label       = "p.format",
    label.y  = ctrl_pos_w)
p_w

## East
# get max value for p position
ctrl_max_e <- df_sweep_all_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_syno, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
  data = df_sweep_all_sel_con %>% filter(pop %in% east_pops),
  x = "pop", 
  y = "de_per_syno",
  color = "region_type", 
  palette = "npg",
  add = "jitter"
) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black") +
  annotate("text", x = 2.5, y=0.8, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 5.5, y=0.8, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Oriental populations", y = NULL, color="Region type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.format",
    label.y = ctrl_pos_e
  )
p_e

p_m <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m

ggsave2("../output/s06.dm_sweep/pop_allele_SNP.region.mutation_burden.svg", p_m, width=10, height=3)
ggsave2("../output/s06.dm_sweep/pop_allele_SNP.region.mutation_burden.pdf", p_m, width=10, height=3)



#########  
## Plot total deleterious numbers
#########


#######################################################################################
#######################################################################################

### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)


# df_sweep <- df_sweep %>%  
#    filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "cauc", "pyra",
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c(
  "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP",
  "White", "pash", "ussu", "betu"
)

df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection")))


### use all_masked regions of the dm
df_sweep_mask <- df_sweep[df_sweep$region_type == "all_masked", ]


# get max value for p position
ctrl_max <- df_sweep_mask %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05


p <- ggboxplot(df_sweep_mask, x = "pop", y = "de_nu",
               color = "type", palette = "aaas") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add p-values
p + stat_compare_means(
  aes(group = type),
  label       = "p.format",
  label.y  = ctrl_pos
) 



#### plot seperately for West and East ####
# keep  All(elle)
cultivar_pops<- c("comm_Dessert", "comm_Perry", "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White")
wild_pops     <- c("pyra", "cauc", "pash", "ussu", "betu")


df_sweep_mask <- df_sweep_mask %>%
  filter(type == "All") %>% 
  mutate(group = ifelse(pop %in% cultivar_pops, "cultivar", "wild"))
  

# get max value for p position
ctrl_max_w <- df_sweep_mask %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# Plot west
p_w <- ggboxplot(
  data = df_sweep_mask %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_nu",
  color = 'group',
  palette = "jco",
  add = "jitter"
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y=10000, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y=10000, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Occidental populations (Ref. P. communis)", y = "No. of Deleterious alleles", color="Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    comparisons = list( c("comm_Dessert", "cauc"),
                       c("comm_Dessert", "pyra"),
                       c("comm_Perry", "cauc"),
                       c("comm_Perry", "pyra")),
    label       = "p.format",
    
    method= 'wilcox.test'
    )
p_w

## East

df_e <- df_sweep_mask %>%
  filter(pop %in% east_pops) %>%
  mutate(group = if_else(pop %in% cultivar_pop_e, "Cultivar", "Wild"))

pval <- wilcox.test(de_nu ~ group, data = df_e)$p.value
y_max <- max(df_e$de_nu, na.rm = TRUE)
x_pos <- mean(c(1, length(unique(df_e$pop)))) # center x position


# get max value for p position
ctrl_max_e <- df_sweep_mask %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_nu, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_mask$pop))) %>%
  pull(max_ctrl) 

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
  data = df_sweep_mask %>% filter(pop %in% east_pops),
  x = "pop", 
  y = "de_nu",
  palette = "jco",
  color = 'group',
  add = "jitter"
) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 2.5, y=10000, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 7, y=10000, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Oriental populations (Ref P. pyrifolia)", y = NULL, color=NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text",
              x    = x_pos,
              y    = y_max * 1.05,
              label= paste0("Wilcox p = ", signif(pval, 3))
           )
p_e
p_m <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m



