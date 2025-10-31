# By Xilong CHEN
# Create date: 2023-09-21
# Contact: chen_xilong@outlook.com

# Modified by Yuqi on 2025-06-11, 06-16, 07-24

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


# Function
save_myplot <- function(p, prefix, w = 180, h = 100, dpi = 600) {
  wi <- w / 25.4
  hi <- h / 25.4
  
  ggsave(paste0(prefix, ".pdf"), p,
         width = wi, height = hi, units = "in")
  
  ggsave(paste0(prefix, ".tif"), p,
         width = wi, height = hi, units = "in",
         dpi = dpi, device = "tiff", compression = "lzw")
  
  ggsave(paste0(prefix, ".png"), p,
         width = wi, height = hi, units = "in", dpi = dpi)
}



###################
## Homo vs He ----
###################

### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)


#df_sweep <- df_sweep %>%  
#    filter(pop  != 'pash' & pop != 'Sand_CN')

# Define population order for Occidental/Eastern

desired_levels <- c(
  "comm_Dessert", "comm_Perry", "pyra", "cauc",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]

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
ctrl_pos_w <- ctrl_max_w * 1.1

# Plot west
p_w <- ggboxplot(
    data = df_sweep_mask %>% filter(pop %in% west_pops),
    x = "pop", 
    y = "de_nu",
    color = "type", 
    palette = "jco",
    add = "jitter", shape = 'type'
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
## per CDS     ----
###################

### read data
rm(list=ls())
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)


# Define population order for Occidental/Eastern
desired_levels <- c(
  "comm_Dessert", "comm_Perry", "pyra", "cauc",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White",
  "pash", "ussu", "betu"
)
west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]
wild_pops     <- c("pyra", "cauc", "pash", "ussu", "betu")


df_sweep <- df_sweep %>%
  mutate(pop = factor(pop, levels = desired_levels)) %>%
  mutate(region_type = factor(region_type, levels = c("all_masked", "control", "selection"))) %>% 
  mutate(crop_type = if_else(pop %in% wild_pops, "Wild", "Cultivar"))
  


df_sweep$de_per_100k_cds <- df_sweep$de_per_cds * 100000


### use all alleles of the dm
df_sweep_all <- df_sweep[df_sweep$type == "All", ]
head(df_sweep_all)

unique(df_sweep_all$region_type)

#### Fig 1 ####
#### use all_masked regions 
df_sweep_all_mask <- df_sweep_all[df_sweep_all$region_type == "all_masked", ]
unique(df_sweep_all_mask$region_type)

my_comparisons <- list( c("comm_Dessert", "pyra"), c("comm_Dessert", "cauc"),
                        c("comm_Dessert", "comm_Perry"), 
                            c("comm_Perry", "cauc"), c("comm_Perry", "pyra")
                        )

#### plot seperately for West and East ####
p_w <- ggboxplot(
  data = df_sweep_all_mask %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "#2400DD",
  fill  = "grey80",  
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 1)
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 10) +
  labs(
    x = "Occidental populations", 
    y = "No. of Deleterious\nalleles/CDS length", 
    color = "Region type",
    fill = "Region type"
  ) +
  annotate("text", x = 1.5, y = 100, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y = 100, label = "wild", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y = 25, label = "", size = 4, fontface = "bold", angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    # aes(group = region_type),
    comparisons = my_comparisons,
    label = "p.format",
   label.y = c(50,60,70,80,90),
   method= 'wilcox.test'
  )

p_w


## East
df_e <- df_sweep_all_mask %>% filter(pop %in% east_pops)

pval <- wilcox.test(de_nu ~ crop_type, data = df_e)$p.value
y_max <- max(df_e$de_nu, na.rm = TRUE)
x_pos <- mean(c(1, length(unique(df_e$pop)))) # center x position


p_e <- ggboxplot(
    data = df_sweep_all_mask %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "#820000", 
    fill = "grey80",
    palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 1)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 3, y=100, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 7, y=100, label = "wild or rootstock", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y = 25, label = "", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Oriental populations", y = NULL, color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  annotate("text",
           x    = x_pos,
           y    = 80,
           label= paste0("Wilcox p=", signif(pval, 3)))
  # ## If add p-values for each group-pair using the following code.
  # stat_compare_means(
  #  #aes(group = region_type),
  #  #comparisons = my_comparisons_e,
  #  label = "p.format"
  # )

p_e

my_comparisons_e <- list(
  c("pash", "Sand_CN"), c("pash", "Sand_CN-SW"), c("pash", "Sand_CN-SE"), c("pash", "pyri_JP"), c("pash", "White"),
  c("ussu", "Sand_CN"), c("ussu", "Sand_CN-SW"), c("ussu", "Sand_CN-SE"), c("ussu", "pyri_JP"), c("ussu", "White"),
  c("betu", "Sand_CN"), c("betu", "Sand_CN-SW"), c("betu", "Sand_CN-SE"), c("betu", "pyri_JP"), c("betu", "White")
)

p_e2 <- ggboxplot(
  data = df_sweep_all_mask %>% filter(pop %in% east_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "#820000", 
  fill = "grey80",
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 1)
)+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 3, y=100, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 7, y=100, label = "wild or rootstock", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y = 25, label = "", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Oriental populations", y = NULL, color="Region type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  # annotate("text",
  #          x    = x_pos,
  #          y    = 80,
  #          label= paste0("Wilcox p=", signif(pval, 3)))
# ## If add p-values for each group-pair using the following code.
stat_compare_means(
 #aes(group = region_type),
 comparisons = my_comparisons_e,
 label = "p.format"
)


p_e2

p_m1 <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_m1

save_myplot(p_m1, "../output/s06.dm_sweep/pop_allele_CDS_1.all_mask.mutation_burden" )



p_m1_2 <- ggarrange(
  p_w, p_e2,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)
p_m1_2

save_myplot(p_m1_2, "../output/s06.dm_sweep/pop_allele_CDS_2_all_p_value.all_mask.mutation_burden" )




#### Fig 2 ####
### use Ho and He alleles of the dm
df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]
unique(df_sweep_all_sel_con$region_type)

my_comparisons <- list( c("all_masked", "control"), c("control", "selection"), c("all_masked", "selection") )

# get max value for p position
ctrl_max <- df_sweep_all_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
  pull(max_ctrl)

# Config color
library(ggsci)
my_cols <- pal_npg()(2)
names(my_cols) <- c("control", "selection")

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

p_w <- ggboxplot(
  data = df_sweep_all_sel_con %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "region_type", 
  fill  = "region_type",  
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 1)
) +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white"
    )
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 100, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 3.5, y = 100, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(
    x = "Occidental populations", 
    y = "No. of Deleterious\nalleles/CDS length", 
    color = "Region type",
    fill = "Region type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    label = "p.signif",
    label.y = ctrl_pos_w
  )
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
    fill = "region_type",
    palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 1)
  )+
  scale_fill_manual(
    values = c(
      "control" = "grey70",
      "selection" = "white")
  ) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 3, y=100, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
  annotate("text", x = 7, y=100, label = "wild", size = 4, fontface = "bold", angle = 0) +
  theme_bw(base_size = 10) +
  labs(x = "Oriental populations", y = NULL, color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.signif",
    method ="wilcox.test",
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




ggsave2("../output/s06.dm_sweep/pop_allele_CDS.selection_and_control.mutation_burden.svg", p_m, width=10, height=3)
ggsave2("../output/s06.dm_sweep/pop_allele_CDS.selection_and_control.all_mask.mutation_burden.pdf", p_m, width=10, height=3)




#### Fig 3 ####
### 0725, 0730
### use Ho and He alleles of the dm
df_sweep_hohe <- df_sweep[df_sweep$type != "All", ]

head(df_sweep_hohe)

df_sweep_hohe_all_mask <- df_sweep_hohe[df_sweep_hohe$region_type == "all_masked", ]

# get max value for p position
ctrl_max <- df_sweep_hohe_all_mask %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_all_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05

p <- ggboxplot(df_sweep_hohe_all_mask, x = "pop", y = "de_per_100k_cds",
               color = "type", palette = "npg") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles/CDS length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

# add p-values
p <- p + stat_compare_means(
    group = type,
    label = "p.format",
    method = "wilcox.test"
    #label.y = ctrl_pos
  )
p

p + stat_compare_means(
  method = "wilcox.test",
  label = "p.format",
  comparisons = list(
    c("control_Ho", "selection_Ho"),   # Only Ho panel
    c("control_He", "selection_He")    # Only He panel
  )
)

p




group_levels <- c(
  "control_Ho", "selection_Ho", 
  "control_He", "selection_He"
)

my_comparisons <- list(
  c("control_Ho", "selection_Ho"),
  c("control_He", "selection_He")
)

df_sweep_hohe_sel_con <- df_sweep_hohe[df_sweep_hohe$region_type != "all_masked", ] %>%
  mutate(group = paste(region_type, type, sep = "_")) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  droplevels()

table(df_sweep_hohe_sel_con$group)



# get max value for p position
ctrl_max <- df_sweep_hohe_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos <- ctrl_max * 1.05

p <- ggboxplot(df_sweep_hohe_sel_con, x = "pop", y = "de_per_100k_cds",
               color = "group", palette = "npg") +
  theme_bw(base_size = 8) +
  labs(x = "Populations", y = "No. of Deleterious alleles/CDS length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(rows = vars(region_type), scales = "free_y", switch = "y") 

# add p-values
p <- p + stat_compare_means(
    comparisons = my_comparisons,
    label = "p.format",
    method = "wilcox.test",
    #label.y = ctrl_pos
  )
p

p + stat_compare_means(
  method = "wilcox.test",
  label = "p.format",
  comparisons = list(
    c("control_Ho", "selection_Ho"),   # Only Ho panel
    c("control_He", "selection_He")    # Only He panel
  )
)

p



#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_hohe_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

npj_colors <- c(
  "control_Ho" = sci_colr(2)[1],
  "selection_Ho" = sci_colr(2)[1],
  "control_He" = sci_colr(2)[2],
  "selection_He" = sci_colr(2)[2]
)

# Plot west
p_w <- ggboxplot(
    data = df_sweep_hohe_sel_con %>% filter(pop %in% west_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "group", 
    fill = "group",
    bxp.errorbar = TRUE,
    palette = "npg",
    add = "jitter",
    add.params = list(alpha = 0.2, size = 1)
  )+
  scale_fill_manual(
    values = c(
      "control_Ho" = "grey80",
      "selection_Ho" = "white",
      "control_He" = "grey80",
      "selection_He" = "white")
  ) +
  scale_color_manual(
    values = njr_olors
  )
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
    annotate("text", x = 1.5, y=80, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 3.5, y=80, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/CDS length", color="Region type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    stat_compare_means(
      aes(group = group),
      label       = "p.format" ) #,
      # label.y  = ctrl_pos_w)
p_w

## East
# get max value for p position
ctrl_max_e <- df_sweep_hohe_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_hohe_sel_con %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "group", 
    palette = "npg",
    add = "jitter"
  ) +
    geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
    annotate("text", x = 3, y=80, label = "cultivar", size = 4, fontface = "bold", angle = 0) +
    annotate("text", x = 7, y=80, label = "wild", size = 4, fontface = "bold", angle = 0) +
    theme_bw(base_size = 10) +
    labs(x = "Oriental populations", y = NULL, color="Region type") 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # stat_compare_means(
  #   #aes(group = region_type),
  #   comparisons = list(c("selection_He", "control_He")),
  #   label = "p.format",
  #   label.y = ctrl_pos_e
  # )

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
  "comm_Dessert", "comm_Perry", "pyra", "cauc",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White",
  "pash", "ussu", "betu"
)

west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]

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
west_pops <- daesired_levels[1:4]
east_pops <- desired_levels[5:12]

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
west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]

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


####NEW 07-25
library(dplyr)
library(ggpubr)
library(ggplot2)

rm(list=ls())

# read data Per SNP 
df_snp <- read.table("../output/s06.dm_sweep/region.SNP.pop_mutation_burden.Rdata.txt", header = T)
df_snp <- df_snp %>%
  mutate(pop = factor(pop, levels = c(
    "comm_Dessert", "comm_Perry", "cauc", "pyra",
    "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
    "pash", "ussu", "betu"
  )),
  region_type = factor(region_type, levels = c("all_masked", "control", "selection"))
  )

west_pops <- c("comm_Dessert", "comm_Perry", "cauc", "pyra")
east_pops <- c("Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White", "pash", "ussu", "betu")

# Only data Ho / He
df_snp_hohe <- df_snp[df_snp$type != "All", ] %>%
  filter(region_type != "all_masked") %>%
  mutate(group = paste(region_type, type, sep = "_"))

# Per SNP west
p_snp_w <- ggboxplot(
  df_snp_hohe %>% filter(pop %in% west_pops),
  x = "pop", y = "de_per_snp", color = "group", palette = "npg", add = "jitter"
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 0.15, label = "cultivar", size = 4, fontface = "bold") +
  annotate("text", x = 3.5, y = 0.15, label = "wild", size = 4, fontface = "bold") +
  labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/SNP") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Per SNP east
p_snp_e <- ggboxplot(
  df_snp_hohe %>% filter(pop %in% east_pops),
  x = "pop", y = "de_per_snp", color = "group", palette = "npg", add = "jitter"
) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 2.5, y = 0.2, label = "cultivar", size = 4, fontface = "bold") +
  annotate("text", x = 7, y = 0.2, label = "wild", size = 4, fontface = "bold") +
  labs(x = "Oriental populations", y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_snp_m <- ggarrange(p_snp_w, p_snp_e, ncol = 2, widths = c(1, 2),
                     common.legend = TRUE, legend = "bottom", align = "hv")


p_snp_m
ggsave2("../output/s06.dm_sweep/pop_HoHe_per_SNP.pdf", p_snp_m, width = 10, height = 3)
ggsave2("../output/s06.dm_sweep/pop_HoHe_per_SNP.svg", p_snp_m, width = 10, height = 3)


# ======================== SYNONYMOUS SNV ===========================

# Read data Synonymous
df_syno <- read.table("../output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt", header = T)
df_syno <- df_syno %>%
  mutate(pop = factor(pop, levels = c(
    "comm_Dessert", "comm_Perry", "cauc", "pyra",
    "Sand_CN", "Sand_CN-SE", "Sand_CN-SW", "pyri_JP", "White",
    "pash", "ussu", "betu"
  )),
  region_type = factor(region_type, levels = c("all_masked", "control", "selection"))
  )

# Only data Ho / He
df_syno_hohe <- df_syno[df_syno$type != "All", ] %>%
  filter(region_type != "all_masked") %>%
  mutate(group = paste(region_type, type, sep = "_"))

# Synonymous west
p_syno_w <- ggboxplot(
  df_syno_hohe %>% filter(pop %in% west_pops),
  x = "pop", y = "de_per_syno", color = "group", palette = "npg", add = "jitter"
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 0.3, label = "cultivar", size = 4, fontface = "bold") +
  annotate("text", x = 3.5, y = 0.3, label = "wild", size = 4, fontface = "bold") +
  labs(x = "Occidental populations", y = "No. of Deleterious\nalleles/synonymous SNPs") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Synonymous east
p_syno_e <- ggboxplot(
  df_syno_hohe %>% filter(pop %in% east_pops),
  x = "pop", y = "de_per_syno", color = "group", palette = "npg", add = "jitter"
) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  annotate("text", x = 2.5, y = 0.4, label = "cultivar", size = 4, fontface = "bold") +
  annotate("text", x = 7, y = 0.4, label = "wild", size = 4, fontface = "bold") +
  labs(x = "Oriental populations", y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_syno_m <- ggarrange(p_syno_w, p_syno_e, ncol = 2, widths = c(1, 2),
                      common.legend = TRUE, legend = "bottom", align = "hv")

p_syno_m
ggsave2("../output/s06.dm_sweep/pop_HoHe_per_synonymous_SNP.pdf", p_syno_m, width = 10, height = 3)
ggsave2("../output/s06.dm_sweep/pop_HoHe_per_synonymous_SNP.svg", p_syno_m, width = 10, height = 3)


