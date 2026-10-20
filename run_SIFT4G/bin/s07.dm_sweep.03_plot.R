# Refered Xilong Chen's script

# Update on 2026-05-28 by Yuqi NIE:
# Using Synonymous sites

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

# stars derived from the p value, rather than hard-coded
signif_star <- function(p) {
  out <- as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                          labels = c("****", "***", "**", "*", "ns")))
  out[is.na(p)] <- ""  # not testable
  out
}

# Per-population paired Wilcoxon signed-rank test between two groups measured on
# the same samples (matched by Sample_ID), then Holm correction within each
# lineage (Occidental / Oriental) separately. Populations with fewer than
# min_n complete pairs are not testable (NA).
paired_pop_tests <- function(dat, group_var, level1, level2,
                             adj_method = "holm", min_n = 3) {
  res <- lapply(unique(as.character(dat$pop)), function(pp) {
    dd <- dat[as.character(dat$pop) == pp, c("Sample_ID", group_var, YVAR)]
    w  <- tidyr::pivot_wider(dd, names_from = all_of(group_var), values_from = all_of(YVAR))
    w  <- w[complete.cases(w[, c(level1, level2)]), ]
    p  <- if (nrow(w) < min_n) NA_real_ else
      tryCatch(suppressWarnings(wilcox.test(w[[level1]], w[[level2]], paired = TRUE)$p.value),
               error = function(e) NA_real_)
    data.frame(pop = pp, p = p, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res$lineage <- ifelse(res$pop %in% west_pops, "Occidental", "Oriental")
  res$p_adj <- NA_real_
  for (lg in unique(res$lineage)) {
    i <- res$lineage == lg
    res$p_adj[i] <- p.adjust(res$p[i], method = adj_method)
  }
  res
}

PAIRED_TEST_NAME <- "two-sided paired Wilcoxon signed-rank test"
adj_note <- function(tests) {
  paste0("Holm, within lineage (", sum(!is.na(tests$p)), " populations tested)")
}

# star labels (from the Holm-adjusted P) at the x/y position of each population
# present in `d`; `ypos` is ordered like levels(d$pop)
star_labels <- function(tests, d, ypos) {
  lv <- levels(d$pop)
  tt <- tests[tests$pop %in% lv, ]
  idx <- match(tt$pop, lv)
  data.frame(xnum = idx, y = ypos[idx], label = signif_star(tt$p_adj))
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else max(x)
}

# y-limit suggestion driven by data. 
# Run once per panel, then paste the value into the YLIM_* constants below.
suggest_ylim <- function(x, probs = c(0.001, 0.999), pad = 0.05) {
  r <- stats::quantile(x, probs, na.rm = TRUE)
  d <- diff(r)
  signif(c(r[1] - d * pad, r[2] + d * pad), 3)
}

# Letter positions clamped to the visible range
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
# P-value bookkeeping, for a unified Supplementary Table
# (Nature Communications requires exact P values to be reported, not just
# significance stars/letters, so every test used to annotate a panel is
# recorded here as it is computed.)
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

# Tidy a pairwise.wilcox.test() p-value matrix into long format, keeping only
# the (non-NA) entries the test actually reports.
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


YVAR      <- "Deleterious_per_synonymous_SNP"
YLAB_MAIN <- "Genetic load\n(per synonymous site)"

SZ_BASE  <- 7          # theme base_size
SZ_ANNOT <- 6 / 2.845  # stars, letters, brackets (pt to size)
SZ_GROUP <- 7 / 2.845  # wild / cultivar labels
JITTER_ALPHA  <- 0.3
JITTER_SIZE   <- 1.5
JITTER_SHPAE  <- 16
JITTER_PARAMS <- list(alpha = JITTER_ALPHA, size = JITTER_SIZE, shape=JITTER_SHPAE)
V_LINE_COLOR <- "grey70"

# Panel a: gap between the significance letter and the "n=" label above it
N_LABEL_OFFSET <- 0.09

# theme_bw without vertical grid lines and without minor horizontal grid lines
theme_pub <- function() {
  theme_bw(base_size = SZ_BASE) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
}

# x positions for two-group paired plots: population index +/- dodge for the
# group, plus a per-sample jitter value.
paired_x <- function(d, group_var, dodge = 0.4, jw = 0.08) {
  ids <- unique(as.character(d$Sample_ID))
  jit <- setNames(runif(length(ids), -jw, jw), ids)
  d %>% mutate(xpos = as.numeric(pop) +
                 (as.numeric(.data[[group_var]]) - 1.5) * dodge +
                 jit[as.character(Sample_ID)])
}

# jittered points, one per sample per group (no connecting line between samples)
paired_layers <- function(d, colour_var) {
  list(
    geom_point(data = d, aes(x = xpos, y = .data[[YVAR]], color = .data[[colour_var]]),
               inherit.aes = FALSE, alpha = JITTER_ALPHA,
               size = JITTER_SIZE, shape = JITTER_SHPAE)
  )
}

# Define population order for Occidental/Eastern
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
YLIM_A <- c(0.15, 0.23); YANNO_A <- 0.210; YBRAK_A <- 0.220   # Fig a
YLIM_B <- c(0.12, 0.45)                                       # Fig b
YLIM_C <- c(0.00, 0.22)                                       # Fig c
YLIM_D <- c(0.00, 0.36)                                       # Fig d
YLIM_S <- c(0.00, 0.35); YANNO_S <- 0.320                     # Fig S

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
df_sweep <- read.table("../output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt", header =T)



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
## Fig a: All (he+ho) alleles ----
###################

# All = Ho + He
df_sweep_all  <- df_sweep[df_sweep$Allelic_state == "All", ]
head(df_sweep_all)

unique(df_sweep_all$region_type)

# all_masked = selection + control
df_sweep_all_mask <- df_sweep_all[df_sweep_all$region_type == "all_masked", ]
unique(df_sweep_all_mask$region_type)

suggest_ylim(df_sweep_all_mask[[YVAR]])
# --- West
df_w <- df_sweep_all_mask %>% filter(pop %in% west_pops)

df_w_test <- df_w %>%
  filter(crop_type %in% c("Cultivar", "Wild")) %>%
  mutate(crop_type = droplevels(factor(crop_type)))
pval_w <- wilcox.test(reformulate("crop_type", YVAR), data = df_w_test)$p.value
message("Fig a West, wild vs cultivar: p = ", signif(pval_w, 3))
record_pval(
  panel = "Fig. 4a", lineage = "Occidental",
  population = "Pooled (pyra+cauc vs comm_Dessert+comm_Perry)",
  comparison = "Wild vs Cultivar",
  group1 = "Wild", group2 = "Cultivar",
  test = "two-sided Wilcoxon rank-sum test",
  p_value = pval_w
)

# pairwise comparisons -> letters
pw        <- pairwise.wilcox.test(df_w[[YVAR]], droplevels(df_w$pop),
                                  p.adjust.method = "holm")
pw_raw    <- pairwise.wilcox.test(df_w[[YVAR]], droplevels(df_w$pop),
                                  p.adjust.method = "none")
M_w       <- tri.to.squ(pw$p.value)
letters_w <- multcompLetters(M_w, threshold = 0.01)$Letters
print(letters_w)

pw_df <- merge(pairwise_matrix_to_df(pw_raw$p.value),
               pairwise_matrix_to_df(pw$p.value),
               by = c("Group1", "Group2"), suffixes = c("_raw", "_holm"))
for (k in seq_len(nrow(pw_df))) {
  record_pval(
    panel = "Fig. 4a", lineage = "Occidental",
    population = paste(pw_df$Group1[k], "vs", pw_df$Group2[k]),
    comparison = "Pairwise population comparison (compact letter display)",
    group1 = pw_df$Group1[k], group2 = pw_df$Group2[k],
    test = "two-sided Wilcoxon rank-sum test",
    p_value = pw_df$P_value_raw[k],
    p_adj = pw_df$P_value_holm[k],
    adj_method = "Holm, within-lineage pairwise comparisons"
  )
}

ypos_w <- letter_y(df_w[[YVAR]], droplevels(df_w$pop), YLIM_A[2])
letters_df_w <- data.frame(
  pop        = names(letters_w),
  letter     = unname(letters_w),
  y_position = as.numeric(ypos_w[names(letters_w)])
)
n_w <- table(droplevels(df_w$pop))
letters_df_w$n_label <- paste0("n=", as.integer(n_w[letters_df_w$pop]))
letters_df_w$n_y     <- letters_df_w$y_position + diff(YLIM_A) * N_LABEL_OFFSET

# West plot
p_w <- ggboxplot(
  data = df_sweep_all_mask %>% filter(pop %in% west_pops),
  x = "pop", 
  y = YVAR,
  color = "pop",
  bxp.errorbar = TRUE,
  add = "jitter",
  add.params = JITTER_PARAMS
) +
  geom_vline(xintercept = VLINE_W, linetype = "dashed", color = V_LINE_COLOR) +
  theme_pub() +
  labs(
    x = NULL,
    y = YLAB_MAIN, 
  ) +
  annotate("text", x = XLAB_W["wild"], y = YANNO_A, label = "Wild", fontface = "bold", size = SZ_GROUP) +
  annotate("text", x = XLAB_W["cultivar"], y = YANNO_A, label = "Cultivar", fontface = "bold", size = SZ_GROUP) +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = w_color)

p_w
 
p_w_letters <- p_w +
  bracket(x1 = XLAB_W["wild"], x2 = XLAB_W["cultivar"],
          y = YBRAK_A, tick = diff(YLIM_A) * 0.02,
          label = signif_star(pval_w), size = SZ_ANNOT) +
  geom_text(data = letters_df_w,
            aes(x = pop, y = y_position, label = letter),
            vjust = 0, size = SZ_ANNOT, inherit.aes = FALSE) +
  geom_text(data = letters_df_w,
            aes(x = pop, y = n_y, label = n_label),
            vjust = 0, size = SZ_ANNOT, colour = V_LINE_COLOR, inherit.aes = FALSE) +
  coord_cartesian(ylim = YLIM_A)

p_w_letters

#######
## --- East
df_e <- df_sweep_all_mask %>% filter(pop %in% east_pops)

df_e_test <- df_e %>%
  filter(crop_type %in% c("Cultivar", "Wild")) %>%
  mutate(crop_type = droplevels(factor(crop_type)))

pval_e <- wilcox.test(reformulate("crop_type", YVAR), data = df_e_test)$p.value
message("Fig a East, wild vs cultivar: p = ", signif(pval_e, 3))
record_pval(
  panel = "Fig. 4a", lineage = "Oriental",
  population = "Pooled (pash+ussu vs Sand_CN+Sand_CN-SW+Sand_CN-SE+pyri_JP+White)",
  comparison = "Wild vs Cultivar",
  group1 = "Wild", group2 = "Cultivar",
  test = "two-sided Wilcoxon rank-sum test",
  p_value = pval_e
)

# significance statistic for east
pe <- pairwise.wilcox.test(df_e[[YVAR]], droplevels(df_e$pop),
                                  p.adjust.method = "holm")
pe_raw <- pairwise.wilcox.test(df_e[[YVAR]], droplevels(df_e$pop),
                                  p.adjust.method = "none")
M_e <- tri.to.squ(pe$p.value)
letters_e <- multcompLetters(M_e, threshold = 0.01)$Letters

pe_df <- merge(pairwise_matrix_to_df(pe_raw$p.value),
               pairwise_matrix_to_df(pe$p.value),
               by = c("Group1", "Group2"), suffixes = c("_raw", "_holm"))
for (k in seq_len(nrow(pe_df))) {
  record_pval(
    panel = "Fig. 4a", lineage = "Oriental",
    population = paste(pe_df$Group1[k], "vs", pe_df$Group2[k]),
    comparison = "Pairwise population comparison (compact letter display)",
    group1 = pe_df$Group1[k], group2 = pe_df$Group2[k],
    test = "two-sided Wilcoxon rank-sum test",
    p_value = pe_df$P_value_raw[k],
    p_adj = pe_df$P_value_holm[k],
    adj_method = "Holm, within-lineage pairwise comparisons"
  )
}

ypos_e <- letter_y(df_e[[YVAR]], droplevels(df_e$pop), YLIM_A[2])
letters_df_e <- data.frame(
  pop        = names(letters_e),
  letter     = unname(letters_e),
  y_position = as.numeric(ypos_e[names(letters_e)])
)
n_e <- table(droplevels(df_e$pop))
letters_df_e$n_label <- paste0("n=", as.integer(n_e[letters_df_e$pop]))
letters_df_e$n_y     <- letters_df_e$y_position + diff(YLIM_A) * N_LABEL_OFFSET

p_e <- ggboxplot(
    data = df_e,
    x = "pop", 
    y = YVAR,
    color = "pop", 
    bxp.errorbar = TRUE,
    add = "jitter",
    add.params = JITTER_PARAMS
  )+
  geom_vline(xintercept = VLINE_E, linetype = "dashed", color = V_LINE_COLOR) +
  annotate("text", x = XLAB_E["wild"], y=YANNO_A, label = "Wild", fontface = "bold", size = SZ_GROUP) +
  annotate("text", x = XLAB_E["rootstock"], y=YANNO_A, label = "Rootstock", fontface = "bold", size = SZ_GROUP) +
  annotate("text", x = XLAB_E["cultivar"], y=YANNO_A, label = "Cultivar", fontface = "bold", size = SZ_GROUP) +
  theme_pub() +
  labs(x = NULL, y = NULL)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none")  +
  scale_color_manual(values = e_color) 

p_e_letters <- p_e +
  bracket(x1 = XLAB_E["wild"], x2 = XLAB_E["cultivar"],
          y = YBRAK_A, tick = diff(YLIM_A) * 0.02,
          label = signif_star(pval_e), size = SZ_ANNOT) +
  geom_text(data = letters_df_e,
            aes(x = pop, y = y_position, label = letter),
            vjust = 0, size = SZ_ANNOT, inherit.aes = FALSE) +
  geom_text(data = letters_df_e,
            aes(x = pop, y = n_y, label = n_label),
            vjust = 0, size = SZ_ANNOT, colour = V_LINE_COLOR, inherit.aes = FALSE) +
  coord_cartesian(ylim = YLIM_A)


p_w_named <- p_w_letters + ggtitle("Occidental pear populations")
p_e_named <- p_e_letters + ggtitle("Oriental pear populations")

p_m1 <- p_w_named + p_e_named +
  plot_layout(ncol = 2, widths = c(1, 2))  &
  theme(plot.title = element_text(hjust = 0.5))  # center

p_m1


###################
## Fig b: Selection vs. Control
###################

### use Ho and He alleles of the dm
# drop all_masked = selection + control
df_sweep_all_sel_con <- df_sweep_all[df_sweep_all$region_type != "all_masked", ] %>%
  droplevels()

suggest_ylim(df_sweep_all_sel_con[[YVAR]]) 

#### plot seperately for West and East ####

# get max value for p position
ctrl_max_w <- df_sweep_all_sel_con %>%
  filter(pop %in% west_pops) %>%
  group_by(pop, .drop = TRUE) %>%
  summarize(max_ctrl = safe_max(.data[[YVAR]]), .groups = "drop") %>%
  arrange(pop) %>%
  pull(max_ctrl)
# * 1.05：
ctrl_pos_w <- ctrl_max_w * 1.05

# paired test per population; Holm within each lineage
tests_b <- paired_pop_tests(df_sweep_all_sel_con, "region_type", "selection", "control")
record_tests <- function(tests, pops, panel, lineage, comparison, group1, group2) {
  tt <- tests[tests$pop %in% pops, ]
  for (k in seq_len(nrow(tt))) {
    record_pval(
      panel = panel, lineage = lineage, population = tt$pop[k],
      comparison = comparison, group1 = group1, group2 = group2,
      test = PAIRED_TEST_NAME,
      p_value = tt$p[k], p_adj = tt$p_adj[k], adj_method = adj_note(tt)
    )
  }
}
record_tests(tests_b, west_pops, "Fig. 4b", "Occidental",
             "Selection vs Control region", "Selection", "Control")

d_2w <- df_sweep_all_sel_con %>% filter(pop %in% west_pops) %>% droplevels() %>%
  paired_x("region_type")

p_2w <- ggboxplot(
  data = d_2w,
  x = "pop",
  y = YVAR,
  color = "region_type",
  palette = "npg",
  bxp.errorbar = TRUE,
  outlier.shape = NA
  ) +
  paired_layers(d_2w, "region_type") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = V_LINE_COLOR) +
  # scale_fill_manual(values = REGION_FILL) +
  # scale_color_manual(values = c(control = "grey40", selection = "black")) +
  theme_pub() +
  labs(x = NULL, 
       y = "Genetic load\n(selection vs. control)",
       color = "Region",
       fill = "Region"
  ) +
  theme(axis.text.x = element_blank()) +
  legend_inside(0.01, 0.98, direction='horizontal') +
  geom_text(data = star_labels(tests_b, d_2w, ctrl_pos_w),
            aes(x = xnum, y = y, label = label),
            size = SZ_ANNOT, inherit.aes = FALSE) +
  coord_cartesian(ylim = YLIM_B)

p_2w

## East
# get max value for p position
ctrl_max_e <- df_sweep_all_sel_con %>%
  filter(pop %in% east_pops) %>%
  group_by(pop, .drop = TRUE) %>%
  summarize(max_ctrl = safe_max(.data[[YVAR]]), .groups = "drop") %>%
  arrange(pop) %>%
  pull(max_ctrl)

# * 1.02
ctrl_pos_e <- pmin(ctrl_max_e * 1.02, YLIM_B[2] * 0.97)

record_tests(tests_b, east_pops, "Fig. 4b", "Oriental",
             "Selection vs Control region", "Selection", "Control")

d_2e <- df_sweep_all_sel_con %>% filter(pop %in% east_pops) %>% droplevels() %>%
  paired_x("region_type")

p_2e <- ggboxplot(
    data = d_2e,
    x = "pop",
    y = YVAR,
    color = "region_type",
    palette = "npg",
    bxp.errorbar = TRUE,
    outlier.shape = NA
  ) +
  paired_layers(d_2e, "region_type") +
  geom_vline(xintercept = VLINE_E, linetype = "dashed", color = V_LINE_COLOR) +
  theme_pub() +
  labs(x = NULL, y = NULL,
      #  color = "Region type",
      #  fill = "Region type"
      ) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none") +
  # legend_inside(0.99, 0.99, direction='horizontal') 
  geom_text(data = star_labels(tests_b, d_2e, ctrl_pos_e),
            aes(x = xnum, y = y, label = label),
            size = SZ_ANNOT, inherit.aes = FALSE) +
  coord_cartesian(ylim = YLIM_B)

p_m2 <- (p_2w + p_2e) + plot_layout(ncol = 2, widths = c(1, 2)) 

p_m2

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
                           panel = "Fig. 4c", lineage = NULL, tests = NULL) {

  d <- dat %>%
    filter(pop %in% pops) %>%
    mutate(pop = factor(pop, levels = pops),
           xnum = as.numeric(pop),
           Allelic_state = factor(Allelic_state, levels = c("Ho", "He")))

  pv <- d %>%
    group_by(pop) %>%
    summarise(
      xnum  = dplyr::first(xnum),
      y_pos = min(safe_max(.data[[YVAR]]) * 1.03, ylim[2] * 0.97),
      .groups = "drop"
    ) %>%
    mutate(p_label = signif_star(tests$p_adj[match(as.character(pop), tests$pop)]))

  record_tests(tests, pops, panel, lineage,
               "Homozygous vs Heterozygous deleterious burden",
               "Homozygous", "Heterozygous")

  ids <- unique(as.character(d$Sample_ID))
  jit <- setNames(runif(length(ids), -0.08, 0.08), ids)
  d_pair <- d %>%
    mutate(xpos = xnum + ifelse(Allelic_state == "Ho", offset, -offset) +
             jit[as.character(Sample_ID)])

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
  # paired lines (Ho - He of the same sample) and jittered points
  paired_layers(d_pair, "pop") +
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
    theme_pub() +
    labs(x = NULL, y = y_lab) +
    theme(axis.text.x = element_blank(),
          axis.text.y = if (show_y_text) element_text() else element_blank())

  p + if (show_legend) legend_inside(0, 0.95) else theme(legend.position = "none")
}

# paired Ho vs He per population; Holm within each lineage
tests_c <- paired_pop_tests(df_sweep_hohe_all_mask, "Allelic_state", "Ho", "He")

p_3w <- zygosity_panel(df_sweep_hohe_all_mask, west_pops, w_color, VLINE_W, YLIM_C,
                      y_lab = "Genetic load\n(by zygosity)",
                      show_legend = FALSE,  # The function show_legend not work yet.
                      panel = "Fig. 4c", lineage = "Occidental", tests = tests_c)
p_3e <- zygosity_panel(df_sweep_hohe_all_mask, east_pops, e_color, VLINE_E, YLIM_C,
                      y_lab = NULL, show_y_text = FALSE,
                      show_legend = FALSE,  # The function show_legend not work yet.
                      panel = "Fig. 4c", lineage = "Oriental", tests = tests_c)

p_m3 <- p_3w + p_3e +
  plot_layout(ncol = 2, widths = c(1, 2))

p_m3



###################
## Fig d: homozygous selection with control
###################
# parepare data
# !All is Ho or He, as All is Ho+He
df_sweep_hohe <- df_sweep[df_sweep$Allelic_state != "All", ]

df_sweep_hohe_sel_con <- df_sweep_hohe %>%
  filter(region_type != "all_masked") %>%
  # rename(zygosity = Allelic_state) %>%   # change column name type to zygosity
  mutate(group = paste(region_type, Allelic_state, sep = "_")) %>%
  droplevels()

df_sweep_hohe_sel_con$group <- factor(df_sweep_hohe_sel_con$group,
                                      levels = unique(df_sweep_hohe_sel_con$group))

df_sweep_ho_sel_con <- df_sweep_hohe_sel_con[df_sweep_hohe_sel_con$Allelic_state == "Ho", ]
df_sweep_he_sel_con <- df_sweep_hohe_sel_con[df_sweep_hohe_sel_con$Allelic_state == "He", ]



#### Fig d main, HO only, and splite West and East
sel_con_panel <- function(dat, pops, vlines, ylim, y_lab = NULL,
                          show_y_text = TRUE, show_legend = FALSE,
                          legend_direction = "vertical",
                          band_labels = NULL, band_y = NULL, mult = 1.05,
                          panel = "Fig. 4d", lineage = NULL, zygosity_label = "Homozygous",
                          tests = NULL) {

  d <- dat %>% filter(pop %in% pops) %>% droplevels()

  lab_pos <- d %>%
    group_by(pop, .drop = TRUE) %>%
    summarize(max_ctrl = safe_max(.data[[YVAR]]), .groups = "drop") %>%
    arrange(pop) %>%
    pull(max_ctrl) * mult
  lab_pos <- pmin(lab_pos, ylim[2] * 0.97)

  record_tests(tests, pops, panel, lineage,
               paste0("Selection vs Control region (", zygosity_label, " sites only)"),
               "Selection", "Control")
  d <- paired_x(d, "region_type")
  p <- ggboxplot(
    data = d, x = "pop", y = YVAR,
    color = "region_type",
    # fill = "region_type",
    palette = "npg",
    bxp.errorbar = TRUE,
    outlier.shape = NA
  )+
  paired_layers(d, "region_type") +
  geom_vline(xintercept = vlines, linetype = "dashed", color = V_LINE_COLOR) +
  scale_x_discrete(labels = SPECIES_LABELS) + 
  # scale_fill_manual(values = REGION_FILL) +
  theme_pub() +
    labs(x = NULL, #"Oriental pear populations",
         y = y_lab,
         color = "Region",
         # fill = "Region"
         ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   colour = unname(pop_color[levels(d$pop)])),
        axis.text.y = if (show_y_text) element_text() else element_blank())+
  geom_text(data = star_labels(tests, d, lab_pos),
            aes(x = xnum, y = y, label = label),
            size = SZ_ANNOT, inherit.aes = FALSE) +
  coord_cartesian(ylim = ylim)

  # p <- p + if (show_legend) legend_inside(0.05, 0.95) else theme(legend.position = "none")
  p <- p + switch(
    as.character(show_legend),
    "TRUE"    = legend_inside(0.01, 0.98, direction = legend_direction),
    "collect" = theme(legend.position = "right"),
    theme(legend.position = "none")
  )

  if (!is.null(band_labels)) {
    for (nm in names(band_labels)) {
      p <- p + annotate("text", x = band_labels[[nm]], y = band_y, label = nm,
                        fontface = "bold", size = SZ_GROUP)
    }
  }
  p
}

# paired selection vs control per population (Ho sites); Holm within each lineage
tests_d <- paired_pop_tests(df_sweep_ho_sel_con, "region_type", "selection", "control")

p_4w <- sel_con_panel(df_sweep_ho_sel_con, west_pops, VLINE_W, YLIM_D,
                     y_lab = "Realized genetic load\n(per synonymous site)",
                     show_legend =TRUE,
                     legend_direction = "horizontal",
                     panel = "Fig. 4d", lineage = "Occidental", zygosity_label = "Homozygous",
                     tests = tests_d)
p_4e <- sel_con_panel(df_sweep_ho_sel_con, east_pops, VLINE_E, YLIM_D,
                     y_lab = NULL, show_y_text = FALSE, mult = 1.05,
                     panel = "Fig. 4d", lineage = "Oriental", zygosity_label = "Homozygous",
                     tests = tests_d)

p_m4 <- (p_4w + p_4e) + plot_layout(ncol = 2, widths = c(1, 2))

p_m4


# p_all4 <- wrap_elements(p_m1) / wrap_elements(p_m2) / wrap_elements(p_m3) / wrap_elements(p_m4) +
#   plot_layout(ncol = 1, heights = c(1, 1, 1, 1.4)) +
#   plot_annotation(tag_levels = 'a') &  # add letter a, b, c...
#   theme(
#     plot.tag = element_text(size = 10, face = "bold")
#   )
# p_all4
# save_myplot(p_all4, "fig4_all4_v3_dm_burden_per_syn_site", w=180, h=170)

p_all4_2 <- p_m1 / p_m2 / p_m3 / p_m4 +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1.4)) +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )
p_all4_2
save_myplot(p_all4_2, "fig4_all4_v3.0.2_dm_burden_per_syn_site", w=180, h=170)


#### Supplementary Figure 1: Potential (He) genetic load, selection vs control

tests_s <- paired_pop_tests(df_sweep_he_sel_con, "region_type", "selection", "control")

p_w <- sel_con_panel(
  df_sweep_he_sel_con, west_pops, VLINE_W, YLIM_S,
  y_lab = "Potential genetic load\n(per synonymous site)",
  show_legend = "collect",
  band_labels = as.list(XLAB_W), band_y = YANNO_S,
  panel = "Supplementary Figure", lineage = "Occidental", zygosity_label = "Heterozygous",
  tests = tests_s
)

p_e <- sel_con_panel(
  df_sweep_he_sel_con, east_pops, VLINE_E, YLIM_S,
  y_lab = NULL, show_y_text = FALSE, mult = 1.02,
  band_labels = as.list(XLAB_E), band_y = YANNO_S,
  panel = "Supplementary Figure", lineage = "Oriental", zygosity_label = "Heterozygous",
  tests = tests_s
)

p_s1 <- (p_w + ggtitle("Occidental pear populations")) +
        (p_e + ggtitle("Oriental pear populations")) +
  plot_layout(ncol = 2, widths = c(1, 2), guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

p_s1

save_myplot(p_s1, "fig_s1_potential_genetic_load_sel_vs_con_per_syn_site", w = 180, h = 70)


#### Supplementary Figure 2: Realized (Ho) over Potential (He) genetic load, not split by selection vs control

YLIM_S2 <- c(0.00, 0.26)
wild_w <- c("pyra", "cauc")
cult_w <- c("comm_Dessert", "comm_Perry")
wild_e <- c("pash", "ussu")
cult_e <- c("Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White")

# single unadjusted Wilcoxon rank-sum test (one comparison per panel), recorded for the Supplementary Table
s2_test <- function(state, state_label, lineage, pops_a, pops_b, group2) {
  d <- df_sweep_hohe_all_mask[df_sweep_hohe_all_mask$Allelic_state == state, ]
  p <- wilcox.test(d[[YVAR]][d$pop %in% pops_a], d[[YVAR]][d$pop %in% pops_b])$p.value
  record_pval(
    panel = "Supplementary Figure", lineage = lineage,
    population = paste0("Pooled (", paste(pops_a, collapse = "+"), " vs ", paste(pops_b, collapse = "+"), ")"),
    comparison = paste0("Wild vs ", group2, " (", state_label, " sites only)"),
    group1 = "Wild", group2 = group2,
    test = "two-sided Wilcoxon rank-sum test",
    p_value = p
  )
  p
}

state_panel <- function(state, pops, vlines, ylim, colorss, grp, brk_to, p_brk,
                        y_lab = NULL, show_y_text = TRUE, show_x_text = FALSE) {
  d <- df_sweep_hohe_all_mask %>%
    filter(Allelic_state == state, pop %in% pops) %>%
    droplevels()
  n_df <- d %>%
    group_by(pop) %>%
    summarise(label = paste0("n=", n()),
              y = max(.data[[YVAR]]) * 1.03 + ylim[2] * 0.02,
              .groups = "drop")
  grp_text <- c(wild = "Wild", rootstock = "Rootstock", cultivar = "Cultivar")[names(grp)]
  ggboxplot(
    data = d, x = "pop", y = YVAR,
    color = "pop",
    bxp.errorbar = TRUE,
    outlier.shape = NA,
    add = "jitter",
    add.params = JITTER_PARAMS
  ) +
    geom_vline(xintercept = vlines, linetype = "dashed", color = V_LINE_COLOR) +
    geom_text(data = n_df, aes(x = pop, y = y, label = label),
              size = SZ_ANNOT, colour = V_LINE_COLOR, vjust = 0, inherit.aes = FALSE) +
    annotate("text", x = unname(grp), y = ylim[2] * 0.9, label = unname(grp_text),
             fontface = "bold", size = SZ_GROUP) +
    bracket(x1 = grp["wild"], x2 = grp[brk_to], y = ylim[2] * 0.95,
            tick = diff(ylim) * 0.02, label = signif_star(p_brk), size = SZ_ANNOT) +
    scale_color_manual(values = colorss) +
    scale_x_discrete(labels = SPECIES_LABELS) +
    theme_pub() +
    labs(x = NULL, y = y_lab) +
    theme(legend.position = "none",
          axis.text.x = if (show_x_text) element_text(angle = 45, hjust = 1, colour = unname(pop_color[levels(d$pop)])) else element_blank(),
          axis.text.y = if (show_y_text) element_text() else element_blank()) +
    coord_cartesian(ylim = ylim)
}

p_ho_w <- s2_test("Ho", "Homozygous", "Occidental", wild_w, cult_w, "Cultivar")
p_ho_e <- s2_test("Ho", "Homozygous", "Oriental", wild_e, cult_e, "Cultivar")
p_he_w <- s2_test("He", "Heterozygous", "Occidental", wild_w, cult_w, "Cultivar")
p_he_e <- s2_test("He", "Heterozygous", "Oriental", wild_e, cult_e, "Cultivar")

s2_hw <- state_panel("Ho", west_pops, VLINE_W, YLIM_S2, w_color, XLAB_W, "cultivar", p_ho_w,
                     y_lab = "Realized genetic load\n(per synonymous site)")
s2_he <- state_panel("Ho", east_pops, VLINE_E, YLIM_S2, e_color, XLAB_E, "cultivar", p_ho_e,
                     show_y_text = FALSE)
s2_pw <- state_panel("He", west_pops, VLINE_W, YLIM_S2, w_color, XLAB_W, "cultivar", p_he_w,
                     y_lab = "Potential genetic load\n(per synonymous site)", show_x_text = TRUE)
s2_pe <- state_panel("He", east_pops, VLINE_E, YLIM_S2, e_color, XLAB_E, "cultivar", p_he_e,
                     show_y_text = FALSE, show_x_text = TRUE)

p_s2 <- (s2_hw + ggtitle("Occidental pear populations")) +
        (s2_he + ggtitle("Oriental pear populations")) +
        s2_pw + s2_pe +
  plot_layout(ncol = 2, widths = c(1, 2)) &
  theme(plot.title = element_text(hjust = 0.5))

p_s2

save_myplot(p_s2, "fig_s2_realized_potential_genetic_load_per_syn_site", w = 180, h = 140)


###################
## Supplementary Table: P values for all panels above
###################

pval_table <- do.call(rbind, pval_records)
pval_table <- pval_table[order(pval_table$Figure_panel, pval_table$Lineage, pval_table$Population), ]
rownames(pval_table) <- NULL

write.table(pval_table, "SupplementaryTable_Fig4_FigS_dm_burden_Pvalues.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
