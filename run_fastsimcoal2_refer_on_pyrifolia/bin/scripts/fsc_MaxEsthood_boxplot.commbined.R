# By Xilong CHEN
# Create date: 2024-11-01
# Contact: chen_xilong@outlook.com

# Update: Yuqi on 2025-10-18

rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(cowplot)   # for theme_half_open() + panel_border()

# -----------------------------
# 1. Global parameters
# -----------------------------
set <- "east"

round_paths <- list(
  round1 = "../../run_fastsimcoal2_refer_on_pyrifolia_round1_4pop/",
  round2 = "../../run_fastsimcoal2_refer_on_pyrifolia_round2b_5pop/",
  roundB1 = "../../run_fastsimcoal2_refer_on_pyrifolia_roundB1b_4pop/",
  roundB2 = "../../run_fastsimcoal2_refer_on_pyrifolia_roundB2_5pop/" #,
  # roundC1 = "../../run_fastsimcoal2_refer_on_pyrifolia_roundC1c1_4pop/"
)

# Output file name for combined figure
combined_outdir <- "./combined_plots/"
dir.create(combined_outdir, showWarnings = FALSE)

# Orange color
highlight_col <- "#E63323"

# -----------------------------
# 2. Define a helper function
# -----------------------------
plot_AIC_boxplot <- function(round_name, round_path, set) {
  
  file_path <- file.path(round_path, "output", "demo_files", 
                         paste0(set, "_summary/MaxEstLhood_AIC_boxplot.txt"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  df <- read.table(file_path, header = TRUE, sep = "\t")
  
  df$divergence <- ifelse(
    nchar(df$Scenario) == 8,
    substr(df$Scenario, nchar(df$Scenario) - 3, nchar(df$Scenario) - 2),
    substr(df$Scenario, nchar(df$Scenario) - 4, nchar(df$Scenario) - 2)
  )
  
  df$gene_flow <- substr(df$Scenario, nchar(df$Scenario) - 1, nchar(df$Scenario)) # last two codes
  
  # Best run per Scenario
  best_runs <- df %>%
    group_by(Scenario) %>%
    slice_min(order_by = AIC, n = 1) %>%
    ungroup()
  
  # Best AIC overall
  best_AIC <- min(best_runs$AIC)
  
  # ΔAIC & Akaike weight
  aic_table <- best_runs %>%
    mutate(
      deltaAIC = AIC - best_AIC,
      weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
    ) %>%
    arrange(AIC)
  
  best_scenario <- aic_table$Scenario[which.min(aic_table$AIC)]
  best_point <- df %>% filter(AIC == min(AIC))
  
  # Build plot
  p <- ggplot(df, aes(x = reorder(Scenario, AIC, FUN = median), y = AIC)) +
    geom_jitter(width = 0.2, alpha = 0.4, color = "grey30") +
    geom_point(data = best_point,
               aes(x = Scenario, y = AIC),
               color = highlight_col, size = 2.5) +
    facet_grid(. ~ gene_flow, scales = "free_x", space = "free_x") +
    theme_half_open(8) + panel_border() +
    labs(
      title = paste0("Round: ", round_name),
      subtitle = paste0("Best scenario: ", best_scenario),
      x = "Scenario",
      y = "AIC"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(face = "bold")
    )
  
  # Save individual plot
  ggsave2(file.path(round_path, "output", "demo_files", 
                    paste0(set, "_summary/MaxEstLhood_boxplot.pdf")), p,
          width = 5, height = 4)
  ggsave2(file.path(round_path, "output", "demo_files", 
                    paste0(set, "_summary/MaxEstLhood_boxplot.png")), p,
          width = 5, height = 4)
  ggsave2(file.path(round_path, "output", "demo_files", 
                    paste0(set, "_summary/MaxEstLhood_boxplot.svg")), p,
          width = 5, height = 4)
  
  return(p)
}

# -----------------------------
# 3. Run the function for all rounds
# -----------------------------
plots <- lapply(names(round_paths), function(rn) {
  plot_AIC_boxplot(rn, round_paths[[rn]], set)
})

# Remove failed rounds (NULLs)
plots <- plots[!sapply(plots, is.null)]

# -----------------------------
# 4. Combine all plots
# -----------------------------
# Arrange plots vertically using patchwork
combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = paste0("AIC distribution across demographic model rounds"),
    theme = theme(plot.title = element_text(face = "bold", size = 8))
  )

# Display combined figure
combined
# Save individual plot
ggsave2(file.path("combined_MaxEstLhood_boxplot.pdf"), combined,
        width = 7.2, height = 9)
ggsave2(file.path("combined_MaxEstLhood_boxplot.png"), combined,
        width = 7.2, height = 9)
ggsave2(file.path("combined_MaxEstLhood_boxplot.svg"), combined,
        width = 7.2, height = 9)






rm(list = ls())
library(ggplot2)
library(cowplot)
library(ggrepel)
library(patchwork)

set='east'

round1_path = '../../run_fastsimcoal2_refer_on_pyrifolia_round1_4pop/'
round2_path = '../../run_fastsimcoal2_refer_on_pyrifolia_round2_5pop/'
roundB1_path = '../../run_fastsimcoal2_refer_on_pyrifolia_roundB1_4pop/'
roundB2_path =  '../../run_fastsimcoal2_refer_on_pyrifolia_roundB2_5pop/'
roundC1_path = '../../run_fastsimcoal2_refer_on_pyrifolia_roundC1b_4pop/'

df <- read.table(paste0(round1_path,"output/demo_files/", set, "_summary/MaxEstLhood_AIC_boxplot.txt"),
    head = T, sep = "\t"
)

df$divgerence <- substr(df$Scenario, 5, 6)
df$gene_flow <- substr(df$Scenario, 7, 8)

library(dplyr)

# Best run for each Scenario
best_runs <- df %>%
  group_by(Scenario) %>%
  slice_min(order_by = AIC, n = 1) %>%  # Run with Minimal AIC
  ungroup()

# Best AIC
best_AIC <- min(best_runs$AIC)

# Cal ΔAIC and AIC weight
aic_table <- best_runs %>%
  mutate(
    deltaAIC = AIC - best_AIC,
    weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
  ) %>%
  arrange(AIC)

aic_table


best_scenario <- aic_table$Scenario[which.min(aic_table$AIC)]
best_point <- df %>% filter(AIC == min(AIC))  # Loacat Best AIC across all runs

# p0 <- ggplot(df, aes(x = Scenario, y = AIC)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#       facet_grid(. ~ gene_flow, scales = "free_x") +
#     theme_half_open(8) +
#     panel_border() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# p0

round1_p0 <- ggplot(df, aes(x = reorder(Scenario, AIC, FUN = median), y = AIC)) +
  #geom_boxplot( outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey30") +
  # Best point
  geom_point(
    data = best_point,
    aes(x = Scenario, y = AIC),
    color = "#E63323", size = 2
  ) +
  facet_grid(. ~ gene_flow, scales = "free_x", space = "free_x") +
  theme_half_open(8) +
  panel_border() +
  labs(
    title = paste0("AIC distribution across models"),
    subtitle = paste0("Best scenario: ", best_scenario),
    x = "Scenario",
    y = "AIC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

round1_p0

# p1 <- ggplot(data = df, aes(x = Scenario, y = MaxEstLhood)) +
#     #geom_boxplot() +
#     facet_grid(. ~ gene_flow, scales = "free_x") +
#     theme_half_open(8) +
#     panel_border() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

ggsave2(paste0(round1_path,"output/demo_files/", set, "_summary/MaxEstLhood_boxplot.pdf"), 
    p0, width = 5, height = 4)
    ggsave2(paste0(round1_path,"output/demo_files/", set, "_summary/MaxEstLhood_boxplot.svg"), 
    p0, width = 5, height = 4)
        ggsave2(paste0(round1_path, "output/demo_files/", set, "_summary/MaxEstLhood_boxplot.png"), 
    p0, width = 5, height = 4)


        
round2_path = '../../run_fastsimcoal2_refer_on_pyrifolia_round2b_5pop/'
        
        df <- read.table(paste0(round2_path,"output/demo_files/", set, "_summary/MaxEstLhood_AIC_boxplot.txt"),
                         head = T, sep = "\t"
        )
        
        df$divgerence <- substr(df$Scenario, 5, 6)
        df$gene_flow <- substr(df$Scenario, 7, 8)
        
        library(dplyr)
        
        # Best run for each Scenario
        best_runs <- df %>%
          group_by(Scenario) %>%
          slice_min(order_by = AIC, n = 1) %>%  # Run with Minimal AIC
          ungroup()
        
        # Best AIC
        best_AIC <- min(best_runs$AIC)
        
        # Cal ΔAIC and AIC weight
        aic_table <- best_runs %>%
          mutate(
            deltaAIC = AIC - best_AIC,
            weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
          ) %>%
          arrange(AIC)
        
        aic_table
        
        
        best_scenario <- aic_table$Scenario[which.min(aic_table$AIC)]
        best_point <- df %>% filter(AIC == min(AIC))  # Loacat Best AIC across all runs
        
        # p0 <- ggplot(df, aes(x = Scenario, y = AIC)) +
        #   geom_boxplot() +
        #   geom_jitter(width = 0.2, alpha = 0.5) +
        #       facet_grid(. ~ gene_flow, scales = "free_x") +
        #     theme_half_open(8) +
        #     panel_border() +
        #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
        # 
        # p0
        
        round2_p0 <- ggplot(df, aes(x = reorder(Scenario, AIC, FUN = median), y = AIC)) +
          #geom_boxplot( outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.4, color = "grey30") +
          # Best point
          geom_point(
            data = best_point,
            aes(x = Scenario, y = AIC),
            color = "#E63323", size = 2
          ) +
          facet_grid(. ~ gene_flow, scales = "free_x", space = "free_x") +
          theme_half_open(8) +
          panel_border() +
          labs(
           # title = paste0("AIC distribution across models"),
            subtitle = paste0("Best scenario: ", best_scenario),
            x = "Scenario",
            y = "AIC"
          ) +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
          )
        
        round2_p0
        
        ggsave2(paste0(round2_path,"output/demo_files/", set, "_summary/MaxEstLhood_boxplot.pdf"), 
                p0, width = 5, height = 4)
        ggsave2(paste0(round2_path,"output/demo_files/", set, "_summary/MaxEstLhood_boxplot.svg"), 
                p0, width = 5, height = 4)
        ggsave2(paste0(round2_path, "output/demo_files/", set, "_summary/MaxEstLhood_boxplot.png"), 
                p0, width = 5, height = 4)

combined <- round1_p0 / round2_p0   # stacked vertically
combined
  
        

# # test G0 vs G1 vs G2 vs G3
# library(ggpubr)
# my_comparisons <- list(c("G1", "G3"))
# ggboxplot(df,
#     x = "gene_flow", y = "MaxEstLhood",
#     color = "gene_flow", palette = "jco"
# ) +
#     stat_compare_means(comparisons = my_comparisons) +
#     stat_compare_means()

# # test 8 models inside G3
# library(multcompView)

# tri.to.squ <- function(x) {
#     rn <- row.names(x)
#     cn <- colnames(x)
#     an <- unique(c(cn, rn))
#     myval <- x[!is.na(x)]
#     mymat <- matrix(1, nrow = length(an), ncol = length(an), dimnames = list(an, an))
#     for (ext in 1:length(cn))
#     {
#         for (int in 1:length(rn))
#         {
#             if (is.na(x[row.names(x) == rn[int], colnames(x) == cn[ext]])) next
#             mymat[row.names(mymat) == rn[int], colnames(mymat) == cn[ext]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
#             mymat[row.names(mymat) == cn[ext], colnames(mymat) == rn[int]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
#         }
#     }
#     return(mymat)
# }

# G3df <- subset(df, df$gene_flow == "G3")

# pp <- pairwise.wilcox.test(G3df$MaxEstLhood, G3df$Scenario, p.adjust.method = "holm")
# mymat <- tri.to.squ(pp$p.value)
# mymat
# myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.01, Letters = letters)
# myletters_df <- data.frame(pop = names(myletters$Letters), letter = myletters$Letters)
# # myletters_df <- myletters_df[-c(6, 7, 8, 9), ]

# p <- ggplot() +
#     geom_jitter(G3df,
#         mapping = aes(x = Scenario, y = MaxEstLhood),
#         width = 0.2, color = "#DDDDDD90", size = 1
#     ) +
#     geom_boxplot(G3df,
#         mapping = aes(x = Scenario, y = MaxEstLhood, fill = Scenario),
#         outlier.shape = NA, alpha = 0.8
#     ) +
#     stat_summary(G3df,
#         mapping = aes(x = Scenario, y = MaxEstLhood, fill = Scenario),
#         fun = mean, geom = "point", shape = 23, size = 1.5
#     ) +
#     theme_half_open(8)
#     theme(legend.position = "none")

# p + geom_text(data = myletters_df, mapping = aes(label = letter, x = pop, y = 2300), colour = "black", size = 3)
