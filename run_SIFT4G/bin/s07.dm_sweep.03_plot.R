# Refered Xilong Chen's script

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(multcompView)
library(patchwork)
library(gghalves) # Required for geom_half_boxplot
library(ggsci)

set.seed(42) # for jitter

# Set working directory to the location of the script and clean up
if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}


# Function save plot, w-max = 180mm, h-max = 170mm
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


# Function transfer triangle df to square df 
tri.to.squ <- function(tri) {
  rows <- rownames(tri)
  cols <- colnames(tri)
  all_groups <- union(rows, cols)

  mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                dimnames = list(all_groups, all_groups))

  for (i in rows) {
    for (j in cols) {
      val <- tri[i, j]
      if (!is.na(val)) {
        mat[i, j] <- val
        mat[j, i] <- val
      }
    }
  }

  return(mat)
}



# 08-22
### read data
df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)

# Define population order for Occidental/Eastern
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
  mutate(pop = factor(pop, levels = desired_levels),
         region_type = factor(region_type, levels = c("all_masked", "control", "selection")),
         crop_type = if_else(pop %in% wild_pops, "Wild", "Cultivar"))

df_sweep$de_per_100k_cds <- df_sweep$de_per_cds * 100000


# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
pop_color <- setNames(color_uniq$Color, color_uniq$Population)

# pops_in_data <- unique(counts_long$Population)
# my_colors <- ori_order_color[pops_in_data]



###################
## Fig a: All (he+ho) alleles ----
###################
# All = Ho + He
df_sweep_all <- df_sweep[df_sweep$type == "All", ]
head(df_sweep_all)

unique(df_sweep_all$region_type)

# all_masked = selection + control
df_sweep_all_mask <- df_sweep_all[df_sweep_all$region_type == "all_masked", ]
unique(df_sweep_all_mask$region_type)

df_w <- df_sweep_all_mask %>% filter(pop %in% west_pops)

# p value for West cultivar and wild
pval <- signif(wilcox.test(de_nu ~ crop_type, data = df_w)$p.value,3)
y_max <- max(df_w$de_per_100k_cds, na.rm = TRUE)
x_pos <- mean(c(1, length(unique(df_w$pop)))) # center x position

pval

p_wc_df <- data.frame(
  group1 = 1.5, # position   
  group2 = 3.5,       # position
  label = pval,
  p_label = "***",     # pval=0.000318
  y.position = y_max*1.3# Position
)


# significance statistic
pw <- pairwise.wilcox.test(df_w$de_per_100k_cds, df_w$pop, p.adjust.method = "holm")
M_w <- tri.to.squ(pw$p.value)
letters_w <- multcompLetters(M_w, threshold = 0.01)$Letters

print(letters_w)
y_pos <- tapply(df_w$de_per_100k_cds, df_w$pop, function(x) max(x, na.rm = TRUE)) * 1.05

letters_df <- data.frame(
  pop = names(letters_w),
  letter = letters_w,
  y_position = y_pos[names(letters_w)])  # make sure the same order of pop in y_pos and letters_w


# West plot
w_color <- pop_color[west_pops]
p_w <- ggboxplot(
  data = df_sweep_all_mask %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "pop",
  # palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 0.6, shape = 16)
) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 7) +
  labs(
    x = NULL, # "Occidental pear populations", 
    y = "No. of deleterious\nalleles / CDS length (100 kb)", 
    color = "Region type",
    fill = "Region type"
  ) +
  annotate("text", x = 1.5, y = 100, label = "cultivar", fontface = "bold", size = 7/2.845,  angle = 0) +
  annotate("text", x = 3.5, y = 100, label = "wild", fontface = "bold", size = 7/2.845,  angle = 0) +
  theme(axis.text.x = element_blank(), # element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_manual(values = w_color)

p_w


p_w_with_pvalue <- p_w + 
  stat_pvalue_manual(
    p_wc_df,
    label = "p_label",
    tip.length = 0.02,
    inherit.aes = FALSE,
    color="grey20",
    size = 6/2.845
  )

p_w_with_pvalue

# add letters
p_w_letters <- p_w_with_pvalue +
  geom_text(data = letters_df,
            aes(x = pop, y = y_position, label = letter),
            vjust = -0.3, size = 6/2.845) +
  coord_cartesian(ylim = c(20, 105))

p_w_letters

#######
## East
df_e <- df_sweep_all_mask %>% filter(pop %in% east_pops)

# significance statistic for east
pe <- pairwise.wilcox.test(df_e$de_per_100k_cds, df_e$pop, p.adjust.method = "holm")
M_e <- tri.to.squ(pe$p.value)
letters_e <- multcompLetters(M_e, threshold = 0.01)$Letters

y_pos <- tapply(df_e$de_per_100k_cds, df_e$pop, function(x) max(x, na.rm = TRUE)) * 1.05

letters_df <- data.frame(
  pop = names(letters_e),
  letter = letters_e,
  y_position = y_pos[names(letters_e)])  # make sure the same order of pop in y_pos and letters_e


# postion
y_max <- max(df_e$de_per_100k_cds, na.rm = TRUE)
x_pos <- mean(c(1, length(unique(df_e$pop)))) # center x position

# p value for East cultivar and wild
pval <- signif(wilcox.test(de_nu ~ crop_type, data = df_e)$p.value,3)
pval


p_wc_df <- data.frame(
  group1 = "Sand_CN-SE", # position   
  group2 = "ussu",       # position
  label = pval,
  p_label = "****",     # pval=1.4e-23
  y.position = y_max*1.2, # Position
  size = 6/2.845
)




e_color <- pop_color[east_pops]
p_e <- ggboxplot(
    data = df_sweep_all_mask %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "pop", 
    # fill = "grey80",
    # palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 0.6, shape=16)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +
  annotate("text", x = 3, y=100, label = "cultivar", fontface = "bold", size = 7/2.845,  angle = 0) +
  annotate("text", x = 7, y=100, label = "wild / rootstock", fontface = "bold", size = 7/2.845,  angle = 0) +
  theme_bw(base_size = 7) +
  labs(x = NULL, y = NULL)+ # "Oriental populations", y = NULL)+ #,  color="Region type") +
    theme(# axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.x = element_blank(),  
          axis.text.y = element_blank(),
          legend.position = "none")  +
  scale_color_manual(values = e_color) 

p_e_with_pvalue <- p_e + 
  stat_pvalue_manual(
    p_wc_df,
    label = "p_label",
    tip.length = 0.02,
    inherit.aes = FALSE,
    color="grey20",
    size = 6/2.845
  )
  
p_e_with_pvalue



p_e_letters <- p_e_with_pvalue+
  geom_text(data = letters_df,
            aes(x = pop, y = y_position, label = letter),
            vjust = -0.3, size = 6/2.845)+
  coord_cartesian(ylim = c(20, 105))

# p_m1 <- ggarrange(
#   p_w_letters, p_e_letters,
#   ncol = 2,
#   nrow = 1,
#   widths = c(1, 2),
#   common.legend = TRUE,
#   legend = "none",
#   align = "hv"
# )

p_w_named <- p_w_letters + ggtitle("Occidental pear populations")
p_e_named <- p_e_letters + ggtitle("Oriental pear populations")

p_m1 <- p_w_named + p_e_named +
  plot_layout(ncol = 2, widths = c(1, 2))  &
  theme(plot.title = element_text(hjust = 0.5))  # center

p_m1

# save_myplot(p_m1, "../output/s06.dm_sweep/pop_allele_CDS_1.all_mask.mutation_burden" )

###################
## Fig b: Selection vs. Control
###################
### use Ho and He alleles of the dm
# drop all_masked = selection + control
df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]
unique(df_sweep_all_sel_con$region_type)

my_comparisons <- list(c("control", "selection") )

## Try facet.by
# df_sweep_all_sel_con <- df_sweep_all_sel_con %>%
#   mutate(
#     group = case_when(
#       pop %in% west_pops ~ "Occidental",
#       pop %in% east_pops ~ "Oriental",
#       TRUE ~ "Error"
#     )
#   )
# ggboxplot(
#      data = df_sweep_all_sel_con,
#      x = "pop", 
#      y = "de_per_100k_cds",
#      facet.by = "group")

# # get max value for p position
# ctrl_max <- df_sweep_all_sel_con %>%
#   group_by(pop) %>%
#   summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
#   arrange(factor(pop, levels = levels(df_sweep_all_sel_con$pop))) %>%
#   pull(max_ctrl)

# Config color
# library(ggsci)
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
  add.params = list(alpha = 0.2, size = 0.6, shape=16)
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white"
    )
  ) +
  theme_bw(base_size = 7) +
  labs(
   x = NULL, # "Occidental populations", 
   y = "Deleterious alleles\n/ CDS (regions)", 
    color = "Region type",
    fill = "Region type"
  ) +
  theme(axis.text.x = element_blank(),
        legend.position = "none"           
        ) +
  stat_compare_means(
    aes(group = region_type),
    label = "p.signif",
    label.y = ctrl_pos_w,
    size = 6/2.845
  )+  coord_cartesian(ylim = c(10,105))

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
    add.params = list(alpha = 0.2, size = 0.6, shape=16)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white")
  ) +
  theme_bw(base_size = 7) +
    labs(x = NULL, # "Oriental populations",
         y = NULL)+ #, color="Region type") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none") +
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.signif",
    method ="wilcox.test",
    label.y = ctrl_pos_e,
    size = 6/2.845
  ) +   coord_cartesian(ylim = c(10,105))

p_m2 <- (p_w + p_e) + plot_layout(ncol = 2, widths = c(1, 2)) & theme(legend.position = 'none')

p_m2

###################
## Fig c: homozygous VS heterozygous states.
###################
### use Ho and He alleles of the dm
# !All is Ho or He, as All is Ho+He
df_sweep_hohe <- df_sweep[df_sweep$type != "All", ]

head(df_sweep_hohe)

df_sweep_hohe_all_mask <- df_sweep_hohe[df_sweep_hohe$region_type == "all_masked", ]


#### plot seperately for West and East ####
### West plot with half boxplot for He
# get max value for p position
ctrl_max_w <- df_sweep_hohe_all_mask %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_all_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# prepare data
df_west <- df_sweep_hohe_all_mask %>% 
  filter(pop %in% west_pops) %>%
  mutate(
    pop = factor(pop, levels = west_pops),
    xnum = as.numeric(pop)
  )

# set size
offset <- 0.2  # Ho/He Horizontal offset winthn a population
w_box <- 0.35   # box wide

# color for pop
w_color <- pop_color[west_pops]

# cal p value
p_values_w <- df_west %>%
  group_by(pop) %>%
  summarise(
    p_val = wilcox.test(de_per_100k_cds[type == "Ho"], 
                        de_per_100k_cds[type == "He"])$p.value,
    xnum = first(xnum),
    y_pos = max(de_per_100k_cds, na.rm = TRUE) + 2,
    .groups = 'drop'
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.0001 ~ "****",
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**", 
      p_val < 0.05 ~ "*",
      p_val < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

# plot
p_w <- ggplot() +
  # Ho: normal boxplot, pop color
geom_boxplot(
  data = filter(df_west, type == "Ho"),
  aes(x = xnum + offset, y = de_per_100k_cds, 
      group = interaction(pop, type), color = pop),
  width = w_box, outlier.shape = NA, outlier.size = 0.8,
  fill = "white"
) +
  stat_boxplot(
    data = filter(df_west, type == "Ho"),
    aes(x = xnum + offset, y = de_per_100k_cds, 
        group = interaction(pop, type), color = pop),
    geom = "errorbar", width = 0.2
  ) +
  # He left: half box, pop color
gghalves::geom_half_boxplot(
  data = filter(df_west, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, 
      group = interaction(pop, "He"), color = pop),
  side = "l", width = w_box, outlier.shape = NA
) +
  # He right: half box, black color
gghalves::geom_half_boxplot(
  data = filter(df_west, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, 
      group = interaction(pop, "He")),
  side = "r", width = w_box, outlier.shape = NA, fill = "black"
) +
  # add jitter for Ho and He
geom_jitter(
  data = filter(df_west, type == "Ho"),
  aes(x = xnum + offset, y = de_per_100k_cds, color = pop),
  width = 0.1, size = 0.6, alpha = 0.2, shape = 16
) +
geom_jitter(
  data = filter(df_west, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, color = pop),
  width = 0.1, size = 0.6, alpha = 0.2, shape = 16
) +
  # add p signif
geom_text(
  data = p_values_w,
  aes(x = xnum, y = y_pos, label = p_label),
  size = 6/2.845, hjust = 0.5
) +
  # split line for wild crop
geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = unique(df_west$xnum), labels = levels(df_west$pop)) +
  scale_color_manual(values = w_color) +
  coord_cartesian(ylim = c(2, 60)) +
  theme_bw(base_size = 7) +
  labs(
    x = NULL,
  y = "Deleterious alleles\n/ CDS (zygosity)",
    # color = "Status"
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )

p_w


## East
# get max value for p position
ctrl_max_e <- df_sweep_hohe_all_mask %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_all_mask$pop))) %>%
  pull(max_ctrl)

# * 1.02：
ctrl_pos_e <- ctrl_max_e * 1.02

# pre data
df_east <- df_sweep_hohe_all_mask %>% 
  filter(pop %in% east_pops) %>%
  mutate(
    pop = factor(pop, levels = east_pops),
    xnum = as.numeric(pop)
  )

# set paramater
offset <- 0.2  # Ho/He horizontal offset
w_box <- 0.35   # box wide

# pop color
e_color <- pop_color[east_pops]

# cal p val
p_values_e <- df_east %>%
  group_by(pop) %>%
  summarise(
    p_val = wilcox.test(de_per_100k_cds[type == "Ho"], 
                        de_per_100k_cds[type == "He"])$p.value,
    xnum = first(xnum),
    y_pos = max(de_per_100k_cds, na.rm = TRUE) + 2,
    .groups = 'drop'
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.0001 ~ "****",
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**", 
      p_val < 0.05 ~ "*",
      p_val < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

# plot
p_e <- ggplot() +
# add jitters for Ho He
geom_jitter(
  data = filter(df_east, type == "Ho"),
  aes(x = xnum + offset, y = de_per_100k_cds, color = pop),
  width = 0.1, size = 0.6, alpha = 0.2, shape = 16
) +
geom_jitter(
  data = filter(df_east, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, color = pop),
  width = 0.1, size = 0.6, alpha = 0.2, shape = 16
) +
  # Ho：noraml box, pop color
geom_boxplot(
  data = filter(df_east, type == "Ho"),
  aes(x = xnum + offset, y = de_per_100k_cds, 
      group = interaction(pop, type), color = pop),
  width = w_box, outlier.shape = NA, outlier.size = 0.8,
  fill = "white"
) +
  stat_boxplot(
    data = filter(df_east, type == "Ho"),
    aes(x = xnum + offset, y = de_per_100k_cds, 
        group = interaction(pop, type), color = pop),
    geom = "errorbar", width = 0.2
  ) +
  # He left box, pop color
gghalves::geom_half_boxplot(
  data = filter(df_east, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, 
      group = interaction(pop, "He"), color = pop),
  side = "l", width = w_box, outlier.shape = NA
) +
  # He right box, pop color
gghalves::geom_half_boxplot(
  data = filter(df_east, type == "He"),
  aes(x = xnum - offset, y = de_per_100k_cds, 
      group = interaction(pop, "He")),
  side = "r", width = w_box, outlier.shape = NA, fill = "black"
) +

  # add p signif
geom_text(
  data = p_values_e,
  aes(x = xnum, y = y_pos, label = p_label),
  size = 6/2.845, hjust = 0.5
) +
  # split line for wild crop
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +

  scale_x_continuous(breaks = unique(df_east$xnum), labels = levels(df_east$pop)) +
  scale_color_manual(values = e_color) +
  coord_cartesian(ylim = c(2, 60)) +
  theme_bw(base_size = 7) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

p_e



# p_m3 <- ggarrange(
#   p_w, p_e,
#   ncol = 2,
#   nrow = 1,
#   widths = c(1, 2),
#   # common.legend = TRUE,
#   # legend = "bottom",
#   align = "hv"
# )
p_m3 <- p_w + p_e +
  plot_layout(ncol = 2, widths = c(1, 2))

p_m3














### 08-26
### Plot not use 
### use Ho and He alleles of the dm
# !All is Ho or He, as All is Ho+He
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

# Preview
p <- ggboxplot(df_sweep_hohe_all_mask, x = "pop", y = "de_per_100k_cds",
               color = "type", palette = "npg") +
  theme_bw(base_size = 7) +
  labs(x = "Populations", y = "Deleterious alleles\n/ CDS (zygosity)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

# add p-values
p_p <- p + stat_compare_means(
    aes(group = .data$type),
    label = "p.signif",
    method = "wilcox.test",
    #label.y=50)
    label.y = ctrl_pos)
p_p

str(p$data)


# Config color
library(ggsci)
my_cols <- pal_npg()(2)
names(my_cols) <- c("control", "selection")

# seperate West and East
### West plot


#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_hohe_all_mask %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_all_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

w_color <- pop_color[west_pops]
p_w <- ggboxplot(
  data = df_sweep_hohe_all_mask %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "type", 
  fill  = "type",  
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 0.6, shape=16)
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "He" = "grey80",
      "Ho" = "white"
    )
  ) +
  theme_bw(base_size = 7) +
  labs(
   x = NULL, # "Occidental populations",
  y = "Deleterious alleles\n/ CDS (zygosity)",
    color = "Status",
    fill = "Status"
  ) +
  theme(axis.text.x = element_blank(),
        legend.position = "none"
        ) +
  stat_compare_means(
    aes(group = type),
    label = "p.signif",
    label.y = ctrl_pos_w,
    size = 6/2.845
  )+  coord_cartesian(ylim = c(4, 60))

p_w


## East
# get max value for p position
ctrl_max_e <- df_sweep_hohe_all_mask %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_all_mask$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02
e_color <- pop_color[east_pops]
p_e <- ggboxplot(
    data = df_sweep_hohe_all_mask %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "type", 
    fill = "type",
    palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 0.6, shape=16)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "He" = "grey70",
      "Ho" = "white")
  ) +
  theme_bw(base_size = 7) +
    labs(x = NULL, # "Oriental populations",
         y = NULL)+ #, color="Region type") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position="none") +
  stat_compare_means(
    aes(group = type),
    #comparisons = list(c("selection", "control")),
    label = "p.signif",
    method ="wilcox.test",
    label.y = ctrl_pos_e,
    size = 6/2.845
  ) +   coord_cartesian(ylim = c(4, 60))

p_e

# p_m3 <- ggarrange(
#   p_w, p_e,
#   ncol = 2,
#   nrow = 1,
#   widths = c(1, 2),
#   common.legend = TRUE,
#   legend = "bottom",
#   align = "hv"
# )






###################
## Fig d: homozygous selection with control
###################

# parepare data
# !All is Ho or He, as All is Ho+He
df_sweep_hohe <- df_sweep[df_sweep$type != "All", ]

df_sweep_hohe_sel_con <- df_sweep_hohe %>%
  filter(region_type != "all_masked") %>%
  rename(zygosity = type) %>%   # change column name type to zygosity
  mutate(group = paste(region_type, zygosity, sep = "_")) %>%
  droplevels()

# for p label position
ctrl_max <- df_sweep_hohe_sel_con %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_hohe_sel_con$pop))) %>%
  pull(max_ctrl)

# ensure the order
group_levels <- unique(df_sweep_hohe_sel_con$group)
df_sweep_hohe_sel_con$group <- factor(df_sweep_hohe_sel_con$group, levels = group_levels)

# Preview 
# plot
p <- ggboxplot(df_sweep_hohe_sel_con,
               x = "pop", y = "de_per_100k_cds",
               color = "group", palette = "npg") +
  theme_bw(base_size = 7) +
  labs(x = "Populations", y = "Deleterious alleles\n/ CDS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(rows = vars(zygosity), scales = "free_y", switch = "y") +
  stat_compare_means(
    aes(group = region_type),        # Selection VS Control
    method = "wilcox.test",
    label = "p.signif",
    label.y = ctrl_max *1.05
  )

p


#### Fig d main, HO only, and splite West and East

# Ho data
head(df_sweep_hohe_sel_con)
df_sweep_ho_sel_con <- df_sweep_hohe_sel_con[df_sweep_hohe_sel_con$zygosity== "Ho",]

head(df_sweep_ho_sel_con)

unique(df_sweep_ho_sel_con$region_type)

# ### use Ho and He alleles of the dm
# # drop all_masked = selection + control
# df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]
# unique(df_sweep_all_sel_con$region_type)

my_comparisons <- list(c("control", "selection") )


# Config color
library(ggsci)
my_cols <- pal_npg()(2)
names(my_cols) <- c("control", "selection")

#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_ho_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_ho_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

p_w <- ggboxplot(
  data = df_sweep_ho_sel_con %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "region_type", 
  fill  = "region_type",  
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 0.6, shape=16)
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white"
    )
  ) +
  theme_bw(base_size = 7) +
  labs(
   x =  NULL, # "Occidental pear populations", 
  y = "Deleterious homozygous\n/ CDS (regions)", 
    color = "Region",
    fill = "Region"
  ) +
  scale_x_discrete(labels = c(
  comm_Dessert = expression(atop(italic("P. communis"), "Dessert")),
  comm_Perry   = expression(atop(italic("P. communis"), "Perry")),
  pyra         = expression(italic("P. pyraster")),
  cauc         = expression(italic("P. caucasica"))
))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.05, 0.95),       
        legend.justification = c(0, 1),        
        legend.background = element_rect(fill = "white", color = "black"),  
        legend.box.background = element_rect(color = "black")               
        ) +
  stat_compare_means(
    aes(group = region_type),
    label = "p.signif",
    label.y = ctrl_pos_w,
    size = 6/2.845
  )+  coord_cartesian(ylim = c(0,80))

p_w


## East
# get max value for p position
ctrl_max_e <- df_sweep_ho_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_ho_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_ho_sel_con %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "region_type", 
    fill = "region_type",
    palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 0.6, shape=16)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white")
  ) +
  theme_bw(base_size = 7) +
    labs(x = NULL, #"Oriental pear populations",
         y = NULL,
         color = "Region",
         fill = "Region"
         ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          legend.position = "none"
          ) +
  scale_x_discrete(labels = c(
  White        = expression(atop(italic("P. pyrifolia"), "White CN")),
  betu         = expression(italic("P. betulifolia")),
  Sand_CN      = expression(atop(italic("P. pyrifolia"), "Sand CN-other")),
  `Sand_CN-SW` = expression(atop(italic("P. pyrifolia"), "Sand CN-SW")),
  `Sand_CN-SE` = expression(atop(italic("P. pyrifolia"), "Sand CN-SE")),
  ussu         = expression(italic("P. ussuriensis")),
  pyri_JP      = expression(atop(italic("P. pyrifolia"), "Nashi JP")),
  pash         = expression(italic("P. pashia"))
))+
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.signif",
    method ="wilcox.test",
    label.y = ctrl_pos_e,
    size = 6/2.845
  ) +   coord_cartesian(ylim = c(0,80))

p_e

# p_m4 <- ggarrange(
#   p_w, p_e,
#   ncol = 2,
#   nrow = 1,
#   widths = c(1, 2),
#   common.legend = TRUE,
#   legend = "bottom",
#   align = "hv"
# )


p_m4 <- (p_w + p_e) + 
  plot_layout(ncol = 2, widths = c(1, 2))


p_m4

# p_all4 <- ggarrange(
#   p_m1, p_m2, p_m3, p_m4,
#   ncol = 1,
#   nrow = 4,
#   # widths = c(1, 2),
#   heights =c(1,1,1,1.8),
#   common.legend = TRUE,
#   legend = "bottom",
#   align = "hv"
# )

p_all4 <- p_m1 / p_m2 / p_m3 / p_m4 +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1)) +
  plot_annotation(tag_levels = 'a') &  # add letter a, b, c...
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )

p_all4


save_myplot(p_all4, "fig4_all4", w=180, h=170)

#### Fig S,  He only, and splite West and East

# He data
head(df_sweep_hohe_sel_con)
df_sweep_he_sel_con <- df_sweep_hohe_sel_con[df_sweep_hohe_sel_con$zygosity== "He",]

head(df_sweep_he_sel_con)

unique(df_sweep_he_sel_con$region_type)

# ### use Ho and He alleles of the dm
# # drop all_masked = selection + control
# df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ]
# unique(df_sweep_all_sel_con$region_type)

my_comparisons <- list(c("control", "selection") )


# Config color
library(ggsci)
my_cols <- pal_npg()(2)
names(my_cols) <- c("control", "selection")

#### plot seperately for West and East ####
# get max value for p position
ctrl_max_w <- df_sweep_he_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_he_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

p_w <- ggboxplot(
  data = df_sweep_he_sel_con %>% filter(pop %in% west_pops),
  x = "pop", 
  y = "de_per_100k_cds",
  color = "region_type", 
  fill  = "region_type",  
  palette = "npg",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = list(alpha = 0.2, size = 0.6, shape=16)
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white"
    )
  ) +
  theme_bw(base_size = 7) +
  labs(
    x =  NULL, # "Occidental pear populations", 
    y = "No. of deleterious heterozygous\nalleles / CDS length (100 kb)", 
    color = "Region",
    fill = "Region"
  ) +
  scale_x_discrete(labels = c(
  comm_Dessert = expression(atop(italic("P. communis"), "Dessert")),
  comm_Perry   = expression(atop(italic("P. communis"), "Perry")),
  pyra         = expression(italic("P. pyraster")),
  cauc         = expression(italic("P. caucasica"))
))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom" # c(0.05, 0.95),       
        # legend.justification = c(0, 1),        
        # legend.background = element_rect(fill = "white", color = "black"),  
        # legend.box.background = element_rect(color = "black")               
        ) +
  annotate("text", x = 1.5, y = 75, label = "cultivar", fontface = "bold", size = 7/2.845,  angle = 0) +
  annotate("text", x = 3.5, y = 75, label = "wild", fontface = "bold", size = 7/2.845,  angle = 0) +
  stat_compare_means(
    aes(group = region_type),
    label = "p.signif",
    label.y = ctrl_pos_w,
    size = 6/2.845
  )+  coord_cartesian(ylim = c(0,85))

p_w


## East
# get max value for p position
ctrl_max_e <- df_sweep_he_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop) %>%
  summarize(max_ctrl = max(de_per_100k_cds, na.rm = TRUE)) %>%
  arrange(factor(pop, levels = levels(df_sweep_he_sel_con$pop))) %>%
  pull(max_ctrl)

# * 1.05：
ctrl_pos_e <- ctrl_max_e * 1.02

p_e <- ggboxplot(
    data = df_sweep_he_sel_con %>% filter(pop %in% east_pops),
    x = "pop", 
    y = "de_per_100k_cds",
    color = "region_type", 
    fill = "region_type",
    palette = "npg",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = list(alpha = 0.2, size = 0.6, shape=16)
  )+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c(
      "control" = "grey80",
      "selection" = "white")
  ) +
  theme_bw(base_size = 7) +
  labs(x = NULL, #"Oriental pear populations",
       y = NULL,
       color = "Region",
       fill = "Region"
      ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "none"
        ) +
  annotate("text", x = 3, y=75, label = "cultivar", fontface = "bold", size = 7/2.845,  angle = 0) +
  annotate("text", x = 7, y=75, label = "wild / rootstock", fontface = "bold", size = 7/2.845,  angle = 0) +
  scale_x_discrete(labels = c(
  White        = expression(atop(italic("P. pyrifolia"), "White CN")),
  betu         = expression(italic("P. betulifolia")),
  Sand_CN      = expression(atop(italic("P. pyrifolia"), "Sand CN-other")),
  `Sand_CN-SW` = expression(atop(italic("P. pyrifolia"), "Sand CN-SW")),
  `Sand_CN-SE` = expression(atop(italic("P. pyrifolia"), "Sand CN-SE")),
  ussu         = expression(italic("P. ussuriensis")),
  pyri_JP      = expression(atop(italic("P. pyrifolia"), "Nashi JP")),
  pash         = expression(italic("P. pashia"))
))+
  stat_compare_means(
    aes(group = region_type),
    #comparisons = list(c("selection", "control")),
    label = "p.signif",
    method ="wilcox.test",
    label.y = ctrl_pos_e,
    size = 6/2.845
  ) +   coord_cartesian(ylim = c(0,85))

p_e



p_w_named <- p_w + ggtitle("Occidental pear populations")
p_e_named <- p_e + ggtitle("Oriental pear populations")

p_m5 <- p_w_named + p_e_named +
  plot_layout(ncol = 2, widths = c(1, 2), guides = "collect") & 
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")  # center

p_m5


save_myplot(p_m5, "fig_s_he_genetic_burden", w=180, h=70)

