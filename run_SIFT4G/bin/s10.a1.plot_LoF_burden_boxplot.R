# Refered Xilong Chen's script

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(multcompView)
library(patchwork)
#  remotes::install_github("erocoar/gghalves")
library(gghalves)  # Required for geom_half_boxplot
library(ggsci)

set.seed(42)  # for jitter

# LoF set from s09: "hc" (high-confidence, main), "consensus" (all putative LoF), "union"
args <- commandArgs(trailingOnly = TRUE)
LOF_SET <- if (length(args) >= 1) args[1] else "hc"

# Set working directory to the location of the script and clean up
if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}

################
# Helper functions
################

# Function save plot, w-max = 180mm, h-max = 170mm
save_myplot <- function(p, prefix, w = 180, h = 100, dpi = 600) {
  wi <- w / 25.4
  hi <- h / 25.4
  
  ggsave(paste0(prefix, ".pdf"), p,
         width = wi, height = hi, units = "in")
  
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

# stars derived from the p value
signif_star <- function(p) {
  out <- as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, Inf),
                          labels = c("****", "***", "**", "*", ".", "ns")))
  out[is.na(p)] <- ""  # not testable
  out
}

# Gguard against groups too small / degenerate for a Wilcoxon test
safe_wilcox_p <- function(x, y, min_n = 3) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < min_n || length(y) < min_n) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(x, y)$p.value), error = function(e) NA_real_)
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else max(x)
}

# y-limit suggestion driven by data
suggest_ylim <- function(x, probs = c(0.001, 0.999), pad = 0.05) {
  r <- stats::quantile(x, probs, na.rm = TRUE)
  d <- diff(r)
  signif(c(r[1] - d * pad, r[2] + d * pad), 3)
}

# Letter positions
letter_y <- function(value, group, ymax, frac = 0.98) {
  pmin(tapply(value, group, safe_max) * 1.02, ymax * frac)
}

# significance bracket
bracket <- function(x1, x2, y, tick, label, size = 6 / 2.845, colour = "black") {
  list(
    annotate("segment", x = x1, xend = x2, y = y, yend = y,
             colour = colour, linewidth = 0.3),
    annotate("segment", x = x1, xend = x1, y = y, yend = y - tick,
             colour = colour, linewidth = 0.3),
    annotate("segment", x = x2, xend = x2, y = y, yend = y - tick,
             colour = colour, linewidth = 0.3),
    annotate("text", x = (x1 + x2) / 2, y = y, label = label,
             vjust = -0.3, size = size, colour = colour)
  )
}

legend_inside <- function(x, y, direction = "vertical") {
  base <- theme(
    legend.justification  = c(0, 1),
    legend.direction      = direction,
    # legend.background     = element_rect(fill = "white", colour = "black"),
    legend.box.background = element_rect(colour = "black")
  )
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    base + theme(legend.position = "inside", legend.position.inside = c(x, y))
  } else {
    base + theme(legend.position = c(x, y))
  }
}

################
# P-value bookkeeping
################

pval_records <- list()

record_pval <- function(panel, lineage, population, comparison,
                         group1, group2,
                         test, p_value, p_adj = NA, adj_method = NA) {
  pval_records[[length(pval_records) + 1]] <<- data.frame(
    Figure_panel = panel,
    Lineage      = lineage,
    Population   = population,
    Comparison   = comparison,
    Group1       = group1,
    Group2       = group2,
    Test         = test,
    P_value      = p_value,
    P_adj        = p_adj,
    Adj_method   = adj_method,
    stringsAsFactors = FALSE
  )
}

# Tidy p-value matrix into long format
pairwise_matrix_to_df <- function(p_mat) {
  rows <- rownames(p_mat)
  cols <- colnames(p_mat)
  out <- do.call(rbind, lapply(rows, function(i) {
    do.call(rbind, lapply(cols, function(j) {
      val <- p_mat[i, j]
      if (is.na(val)) return(NULL)
      data.frame(Group1 = i, Group2 = j, P_value = val, stringsAsFactors = FALSE)
    }))
  }))
  out
}


YVAR      <- "LoF_per_synonymous_site"
YLAB_MAIN <- "LoF burden\n(per synonymous site)"

SZ_BASE  <- 7          # theme base_size
SZ_ANNOT <- 6 / 2.845  # stars, letters, brackets (pt to size)
SZ_GROUP <- 7 / 2.845  # wild / cultivar labels
JITTER_ALPHA  <- 0.3
JITTER_SIZE   <- 1.5
JITTER_SHPAE  <- 16
JITTER_PARAMS <- list(alpha = JITTER_ALPHA, size = JITTER_SIZE, shape=JITTER_SHPAE)
V_LINE_COLOR <- "grey70"

# gap between significance letter and "n="
N_LABEL_OFFSET <- 0.09

# Population order for Occidental/Oriental
desired_levels <- c(
  "pyra", "cauc",
  "comm_Dessert", "comm_Perry",
  "pash", "ussu", "betu",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White"
)

west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]
wild_pops <- c("pyra", "cauc", "pash", "ussu")
rootstock_pops <- c("betu")
crop_pops      <- c("comm_Dessert", "comm_Perry",
                    "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White")

# Panel lines and labels position
VLINE_W <- 2.5
XLAB_W  <- c(wild = 1.5, cultivar = 3.5)

# East x order: pash(1) ussu(2) | betu(3) | Sand_CN(4) ...
VLINE_E <- c(2.5, 3.5)
XLAB_E  <- c(wild = 1.5, rootstock = 3, cultivar = 6)

# --- Y limits
YLIM_A <- c(0.014, 0.034); YANNO_A <- 0.0302; YBRAK_A <- 0.0318   # Fig a (total)
YLIM_C <- c(0.000, 0.018)                                         # Fig c (by zygosity), not plotted
YLIM_HO <- c(0.000, 0.0246); YANNO_HO <- 0.0194; YBRAK_HO <- 0.0220  # Fig c (realized, Ho)
YLIM_HE <- c(0.0078, 0.0200); YANNO_HE <- 0.0176; YBRAK_HE <- 0.0188 # Fig e (potential, He)

# HC set
if (LOF_SET == "hc") {
  YLIM_A <- c(0.0085, 0.0172); YANNO_A <- 0.0156; YBRAK_A <- 0.0164   # Fig a (total), room for "n=" labels
  YLIM_C <- c(0.000, 0.0095)                                          # Fig c (by zygosity), not plotted
  YLIM_HO <- c(0.000, 0.0110); YANNO_HO <- 0.0087; YBRAK_HO <- 0.0099   # Fig c (realized, Ho)
  YLIM_HE <- c(0.0051, 0.0105); YANNO_HE <- 0.0094; YBRAK_HE <- 0.0099  # Fig e (potential, He)
}

# Italic species on x-axis
SPECIES_LABELS <- c(
  pyra         = expression(italic("P. pyraster")),
  cauc         = expression(italic("P. caucasica")),
  comm_Dessert = expression(atop(italic("P. communis"), "Dessert")),
  comm_Perry   = expression(atop(italic("P. communis"), "Perry")),
  pash         = expression(italic("P. pashia")),
  ussu         = expression(italic("P. ussuriensis")),
  betu         = expression(italic("P. betulifolia")),
  Sand_CN      = expression(atop(italic("P. pyrifolia"), "Sand CN-other")),
  `Sand_CN-SW` = expression(atop(italic("P. pyrifolia"), "Sand CN-SW")),
  `Sand_CN-SE` = expression(atop(italic("P. pyrifolia"), "Sand CN-SE")),
  pyri_JP      = expression(atop(italic("P. pyrifolia"), "Nashi JP")),
  White        = expression(atop(italic("P. pyrifolia"), "White CN"))
)

REGION_FILL <- c(control = "grey80", selection = "white")

### read data
# df_sweep <- read.table("../output/s06.dm_sweep/region.CDS.pop_mutation_burden.Rdata.txt", header =T)
df_sweep <- read.table(paste0("../output/s09.LoF_burden/", LOF_SET, "/region.synonymous.pop_mutation_burden.Rdata.txt"), header =T)



df_sweep <- df_sweep %>%
  mutate(pop = factor(Population, levels = desired_levels),
         region_type = factor(Region_type, levels = c("all_masked", "control", "selection")),
         # crop_type = if_else(pop %in% wild_pops, "Wild", "Cultivar"))
         crop_type = case_when(pop %in% wild_pops ~ "Wild",
                               pop %in% rootstock_pops ~ "Rootstock",
                               pop %in% crop_pops ~ "Cultivar",
                               TRUE ~ "other"))


# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config %>% distinct(Population, Color)
pop_color <- setNames(color_uniq$Color, color_uniq$Population)

w_color <- pop_color[west_pops]
e_color <- pop_color[east_pops]


###################
## Fig a,b: All (he+ho),
#      c,d: Realized (Ho),
#      e,f: Potential (He)
###################

GROUP_LABELS <- c(wild = "Wild", rootstock = "Rootstock", cultivar = "Cultivar")

burden_panel <- function(dat, pops, colorss, vlines, xlab_pos, ylim, yanno, ybrak,
                         y_lab = NULL, show_y_text = TRUE, show_n = FALSE,
                         panel = NULL, lineage = NULL, pooled_label = NULL) {

  d <- dat %>% filter(pop %in% pops)

  # wild vs cultivar
  d_test <- d %>%
    filter(crop_type %in% c("Cultivar", "Wild")) %>%
    mutate(crop_type = droplevels(factor(crop_type)))
  pval <- wilcox.test(reformulate("crop_type", YVAR), data = d_test)$p.value
  message(panel, " ", lineage, ", wild vs cultivar: p = ", signif(pval, 3))
  record_pval(
    panel = panel, lineage = lineage, population = pooled_label,
    comparison = "Wild vs Cultivar",
    group1 = "Wild", group2 = "Cultivar",
    test = "two-sided Wilcoxon rank-sum test",
    p_value = pval
  )

  # pairwise comparisons -> letters
  pw     <- pairwise.wilcox.test(d[[YVAR]], droplevels(d$pop), p.adjust.method = "holm")
  pw_raw <- pairwise.wilcox.test(d[[YVAR]], droplevels(d$pop), p.adjust.method = "none")
  lett   <- multcompLetters(tri.to.squ(pw$p.value), threshold = 0.01)$Letters
  print(lett)

  pw_df <- merge(pairwise_matrix_to_df(pw_raw$p.value),
                 pairwise_matrix_to_df(pw$p.value),
                 by = c("Group1", "Group2"), suffixes = c("_raw", "_holm"))
  for (k in seq_len(nrow(pw_df))) {
    record_pval(
      panel = panel, lineage = lineage,
      population = paste(pw_df$Group1[k], "vs", pw_df$Group2[k]),
      comparison = "Pairwise population comparison (compact letter display)",
      group1 = pw_df$Group1[k], group2 = pw_df$Group2[k],
      test = "two-sided Wilcoxon rank-sum test",
      p_value = pw_df$P_value_raw[k],
      p_adj = pw_df$P_value_holm[k],
      adj_method = "Holm, within-lineage pairwise comparisons"
    )
  }

  ypos <- letter_y(d[[YVAR]], droplevels(d$pop), ylim[2])
  letters_df <- data.frame(
    pop        = names(lett),
    letter     = unname(lett),
    y_position = as.numeric(ypos[names(lett)])
  )
  n_pop <- table(droplevels(d$pop))
  letters_df$n_label <- paste0("n=", as.integer(n_pop[letters_df$pop]))
  letters_df$n_y     <- letters_df$y_position + diff(ylim) * N_LABEL_OFFSET

  p <- ggboxplot(
    data = d,
    x = "pop",
    y = YVAR,
    color = "pop",
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = JITTER_PARAMS
  ) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = V_LINE_COLOR) +
    theme_bw(base_size = SZ_BASE) +
    labs(x = NULL, y = y_lab) +
    theme(axis.text.x = element_blank(),
          axis.text.y = if (show_y_text) element_text() else element_blank(),
          legend.position = "none") +
    scale_color_manual(values = colorss)

  for (nm in names(xlab_pos)) {
    p <- p + annotate("text", x = xlab_pos[[nm]], y = yanno, label = GROUP_LABELS[[nm]],
                      fontface = "bold", size = SZ_GROUP)
  }

  p <- p +
    bracket(x1 = xlab_pos["wild"], x2 = xlab_pos["cultivar"],
            y = ybrak, tick = diff(ylim) * 0.02,
            label = signif_star(pval), size = SZ_ANNOT) +
    geom_text(data = letters_df,
              aes(x = pop, y = y_position, label = letter),
              vjust = 0, size = SZ_ANNOT, inherit.aes = FALSE)

  if (show_n) {
    p <- p +
      geom_text(data = letters_df,
                aes(x = pop, y = n_y, label = n_label),
                vjust = 0, size = SZ_ANNOT, colour = V_LINE_COLOR, inherit.aes = FALSE)
  }

  p + coord_cartesian(ylim = ylim)
}

POOLED_W <- "Pooled (pyra+cauc vs comm_Dessert+comm_Perry)"
POOLED_E <- "Pooled (pash+ussu vs Sand_CN+Sand_CN-SW+Sand_CN-SE+pyri_JP+White)"

# Occidental + Oriental panels of one allelic state ("All", "Ho" or "He")
burden_row <- function(state, ylim, yanno, ybrak, y_lab, panel, show_n = FALSE) {
  dat <- df_sweep %>% filter(Allelic_state == state, region_type == "all_masked")
  print(suggest_ylim(dat[[YVAR]]))

  p_w <- burden_panel(dat, west_pops, w_color, VLINE_W, XLAB_W, ylim, yanno, ybrak,
                      y_lab = y_lab, show_n = show_n, panel = panel, lineage = "Occidental",
                      pooled_label = POOLED_W)
  p_e <- burden_panel(dat, east_pops, e_color, VLINE_E, XLAB_E, ylim, yanno, ybrak,
                      y_lab = NULL, show_y_text = FALSE, show_n = show_n,
                      panel = panel, lineage = "Oriental", pooled_label = POOLED_E)
  list(w = p_w, e = p_e)
}

# a,b: total LoF burden
row_all <- burden_row("All", YLIM_A, YANNO_A, YBRAK_A,
                      y_lab = YLAB_MAIN, panel = "Fig. S_LoF a-b", show_n = TRUE)
# c,d: realized LoF burden (homozygous)
row_ho  <- burden_row("Ho", YLIM_HO, YANNO_HO, YBRAK_HO,
                      y_lab = "Realized LoF burden\n(per synonymous site)", panel = "Fig. S_LoF c-d")
# e,f: potential LoF burden (heterozygous)
row_he  <- burden_row("He", YLIM_HE, YANNO_HE, YBRAK_HE,
                      y_lab = "Potential LoF burden\n(per synonymous site)", panel = "Fig. S_LoF e-f")

p_m1 <- (row_all$w + ggtitle("Occidental pear populations")) +
        (row_all$e + ggtitle("Oriental pear populations")) +
  plot_layout(ncol = 2, widths = c(1, 2)) &
  theme(plot.title = element_text(hjust = 0.5))  # center

p_m1

###################
## Fig c: homozygous VS heterozygous states.
###################

### use Ho and He alleles of the dm
# !All is Ho or He, as All is Ho+He
df_sweep_hohe <- df_sweep[df_sweep$Allelic_state != "All", ]

head(df_sweep_hohe)

df_sweep_hohe_all_mask <- df_sweep_hohe[df_sweep_hohe$region_type == "all_masked", ]

suggest_ylim(df_sweep_hohe_all_mask[[YVAR]])

offset <- 0.2   # Ho/He horizontal offset within a population
w_box  <- 0.35  # box width

zygosity_panel <- function(dat, pops, colorss, vlines, ylim, y_lab = NULL,
                           show_y_text = TRUE, show_legend = FALSE,
                           panel = "Fig. S_LoF c-d", lineage = NULL) {

  d <- dat %>%
    filter(pop %in% pops) %>%
    mutate(pop = factor(pop, levels = pops),
           xnum = as.numeric(pop),
           Allelic_state = factor(Allelic_state, levels = c("Ho", "He")))

  pv <- d %>%
    group_by(pop) %>%
    summarise(
      p_val = safe_wilcox_p(.data[[YVAR]][Allelic_state == "Ho"],
                            .data[[YVAR]][Allelic_state == "He"]),
      xnum  = dplyr::first(xnum),
      y_pos = min(safe_max(.data[[YVAR]]) * 1.03, ylim[2] * 0.97),
      .groups = "drop"
    ) %>%
    mutate(
          #  xum = as.numeric(pop),
           p_label = signif_star(p_val))

  for (k in seq_len(nrow(pv))) {
    record_pval(
      panel = panel, lineage = lineage, population = as.character(pv$pop[k]),
      comparison = "Homozygous vs Heterozygous LoF burden",
      group1 = "Homozygous (Ho)", group2 = "Heterozygous (He)",
      test = "two-sided Wilcoxon rank-sum test",
      p_value = pv$p_val[k]
    )
  }

  p <- ggplot() +
  # Ho: normal full boxplot, pop color
    geom_boxplot(
      data = filter(d, Allelic_state == "Ho"),
      aes(x = xnum + offset, y = .data[[YVAR]], 
        group = interaction(pop, Allelic_state), 
        color = pop, fill = Allelic_state),
    width = w_box, outlier.shape = NA
    ) +
  stat_boxplot(
    data = filter(d, Allelic_state == "Ho"),
    aes(x = xnum + offset, y = .data[[YVAR]], 
        group = interaction(pop, Allelic_state), color = pop),
    geom = "errorbar", width = 0.2
  ) +
  # He left: half box filled pop color
gghalves::geom_half_boxplot(
  data = filter(d, Allelic_state == "He"),
  aes(x = xnum - offset, y = .data[[YVAR]], 
      group = interaction(pop, Allelic_state), color = pop),
  side = "l", width = w_box, outlier.shape = NA
) +
  # He right: half box, black color
gghalves::geom_half_boxplot(
  data = filter(d, Allelic_state == "He"),
  aes(x = xnum - offset, y = .data[[YVAR]], 
      group = interaction(pop, Allelic_state)),
  side = "r", width = w_box, outlier.shape = NA, fill = "black"
) +
  # add jitter for Ho and He
geom_jitter(
  data = filter(d, Allelic_state == "Ho"),
  aes(x = xnum + offset, y = .data[[YVAR]], color = pop),
  width = 0.1, size = JITTER_SIZE, alpha = JITTER_ALPHA, shape = JITTER_SHPAE
) +
geom_jitter(
  data = filter(d, Allelic_state == "He"),
  aes(x = xnum - offset, y = .data[[YVAR]], color = pop),
  width = 0.1, size = JITTER_SIZE, alpha = JITTER_ALPHA, shape = JITTER_SHPAE
) +
  # add p signif
geom_text(
  data = pv,
  aes(x = xnum, y = y_pos, label = p_label),
  size = SZ_ANNOT, hjust = 0.5, inherit.aes =FALSE) +
  # split line for wild crop
  geom_vline(xintercept = vlines, linetype = "dashed", color = V_LINE_COLOR) +
  scale_x_continuous(breaks = seq_along(levels(d$pop)), labels = levels(d$pop)) +
  # scale_x_continuous(breaks = unique(d$xnum), labels = levels(d$pop)) +
  scale_color_manual(values = colorss) +
  scale_fill_manual(
      values = c(Ho = "white", He = "black"),
      breaks = c("Ho", "He"),
      labels = c("Homozygous", "Heterozygous"),
      name   = NULL
    ) +
    guides(color = "none",
           fill  = guide_legend(override.aes = list(colour = "grey30"))) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = SZ_BASE) +
    labs(x = NULL, y = y_lab) +
    theme(axis.text.x = element_blank(),
          axis.text.y = if (show_y_text) element_text() else element_blank())

  p + if (show_legend) legend_inside(0, 0.95) else theme(legend.position = "none")
}

###################
## Combine: genome-wide LoF burden only
###################

p_m2 <- row_ho$w + row_ho$e +
  plot_layout(ncol = 2, widths = c(1, 2))

# e,f add species labels on the x axis
p_m3_x <- (row_he$w +
             scale_x_discrete(labels = SPECIES_LABELS[west_pops]) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))) +
          (row_he$e +
             scale_x_discrete(labels = SPECIES_LABELS[east_pops]) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))) +
  plot_layout(ncol = 2, widths = c(1, 2))

p_all2 <- p_m1 / p_m2 / p_m3_x +
  plot_layout(ncol = 1, heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )
p_all2
if (LOF_SET != "consensus") save_myplot(p_all2, paste0("figS_LoF_", LOF_SET, "_genome_wide_burden_per_syn_site"), w = 180, h = 160)


###################
## Supplementary Table: P values for all panels above
###################

pval_table <- do.call(rbind, pval_records)
pval_table <- pval_table[order(pval_table$Figure_panel, pval_table$Lineage, pval_table$Population), ]
rownames(pval_table) <- NULL

write.table(pval_table, paste0("SupplementaryTable_FigS_LoF_", LOF_SET, "_burden_Pvalues.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
