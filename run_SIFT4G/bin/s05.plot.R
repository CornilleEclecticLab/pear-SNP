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
    header = TRUE,
    col.names = c("ind", "pop", "type", "de_nu", "tol_nu", "syn_nu")
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
    mutate(
        Type = type)




# Plot Allele ------------------------------------------------------------
dir.create("../plot", showWarnings = FALSE)
de_all <- de_pop[de_pop$type == "Allele", ]

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
de_hohe <- de_pop[de_pop$Type != "Allele", ]

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
de_he <- de_pop[de_pop$type == "Heterozygous", ]

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
de_ho <- de_pop[de_pop$type == "Homozygous", ]

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