library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(multcompView)
library(ggpubr)

# Set working directory to the location of the script and clean up
if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list = ls())
}

# Load data
de_pop <- read.table(
    "../output/s04.analysis_mutation_burden/pop_mutation_burden.txt",
    header = TRUE
)

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

de_pop$pop <- factor(de_pop$pop, levels = desired_levels)

# Calculate deleterious/tolerated ratio
de_pop <- de_pop %>%
    mutate(del_tol = de_nu / tol_nu) 

de_pop <- de_pop %>% 
    rename(
        Type = type)




# Plot Allele ------------------------------------------------------------
de_all <- de_pop[de_pop$type == "All", ]

p <- ggplot() +
    geom_jitter(
        de_all,
        mapping = aes(x = pop, y = del_tol),
        width = 0.2,
        color = "#DDDDDD90",
        size = 1
    ) +
    geom_boxplot(
        de_all,
        mapping = aes(x = pop, y = del_tol, fill = pop),
        outlier.shape = NA
    ) +
    stat_summary(
        de_all,
        mapping = aes(x = pop, y = del_tol, fill = pop),
        fun = mean,
        geom = "point",
        shape = 23,
        size = 1.5
    ) +
    theme_bw(8) +
    theme(legend.position = "none")

p

# Perform pairwise Wilcoxon test
pp <- pairwise.wilcox.test(de_all$de_nu, de_all$pop, p.adjust.method = "holm")
mymat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.01, Letters = letters)

myletters_df <- data.frame(
    pop = names(myletters$Letters),
    letter = myletters$Letters
)

p + geom_text(
    data = myletters_df,
    mapping = aes(label = letter, x = pop, y = 0.3),
    colour = "black",
    size = 3
)

ggsave2(
    "../plot/pop_all_deleterious_mutation_box_plot.Mar24.svg",
    width = 2.25,
    height = 2.5
)



# Homo Hete ---------------------------------------------------------------
de_hohe <- de_pop[de_pop$Type != "All", ]

p_w <- ggboxplot(
  data = de_hohe %>% filter(pop %in% west_pops),
  x = "pop",
  y = "de_nu",
  color = "Type",
  palette = "npg"
) +
  labs(
    x = "Occidental Populations",
    y = "Deleterious number"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  stat_compare_means(aes(group = Type), label = "p.signif"
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black")



p_e <- ggboxplot(
  data = de_hohe %>% filter(pop %in% east_pops),
  x = "pop",
  y = "de_nu",
  color = "Type",
  palette = "npg"
) +
  labs(
    x = "Oriental Populations",
    y = element_blank()) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Type))+#label = "p.signif") +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")

p_c <- ggarrange(
  p_w, p_e,
  ncol = 2,
  nrow = 1,
  widths = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

p_c


pdfname <- "../plot/pop_HevsHo_deleterious_mutation_box_plot.pdf"
ggsave2(pdfname, device = "pdf", width = 5.5, height = 4)

# He ----------------------------------------------------------------------
de_he <- de_pop[de_pop$type == "He", ]

p <- ggplot(de_he, aes(x = pop, y = de_nu, fill = pop)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3) +
    theme_classic() +
    theme(legend.position = "none")

pp <- pairwise.wilcox.test(de_he$de_nu, de_he$pop, p.adjust.method = "holm")
mymat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.01, Letters = letters)

myletters_df <- data.frame(
    pop = names(myletters$Letters),
    letter = myletters$Letters
)

p + geom_text(
    data = myletters_df,
    aes(label = letter, y = 2000),
    colour = "black",
    size = 4
)

ggsave2("../plot/pop_He_deleterious_mutation_box_plot.pdf", width = 4, height = 3)

# Ho ----------------------------------------------------------------------
de_ho <- de_pop[de_pop$type == "Ho", ]

p <- ggplot(de_ho, aes(x = pop, y = de_nu, fill = pop)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3) +
    theme_classic() +
    theme(legend.position = "none")

pp <- pairwise.wilcox.test(de_ho$de_nu, de_ho$pop, p.adjust.method = "holm")
mymat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.01, Letters = letters)

myletters_df <- data.frame(
    pop = names(myletters$Letters),
    letter = myletters$Letters
)

p + geom_text(
    data = myletters_df,
    aes(label = letter, y = 800),
    colour = "black",
    size = 4
)

ggsave2("../plot/pop_Ho_deleterious_mutation_box_plot.pdf", width = 4, height = 3)