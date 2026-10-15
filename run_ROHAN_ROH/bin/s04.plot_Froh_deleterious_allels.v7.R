library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(ggpubr)
library(patchwork)

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

  ggsave(paste0(prefix, ".png"), p,
    width = wi, height = hi, units = "in", dpi = dpi)
}

# Loading color palette for each population
color_config <- read.table("../input/population_color.tsv", header = TRUE)
color_uniq <- color_config |> distinct(Population, Color)
pop_color <- setNames(color_uniq$Color, color_uniq$Population)

type_mapping <- c(
  "ussu" = "Wild", "pash" = "Wild", "cauc" = "Wild", "pyra" = "Wild",
  "betu" = "Rootstock",
  "comm_Dessert" = "Cultivar", "comm_Perry" = "Cultivar", "pyri_JP" = "Cultivar",
  "Sand_CN" = "Cultivar", "Sand_CN-SE" = "Cultivar", "Sand_CN-SW" = "Cultivar", "White" = "Cultivar"
)

# 17: triangle,    7:box+x;         4:x 
type_shape <- c("Wild" = 17, "Rootstock" = 7, "Cultivar" = 4)

# Load FROH data (Occidental & Oriental) and SIFT data
file_froh_w <- "../../run_ROHAN_ROH_refer_on_communis/output/FrohMatrix.tsv"
file_froh_e <- "../../run_ROHAN_ROH/output/FrohMatrix.tsv"
file_sift   <- "../../run_SIFT4G/output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt"

df_froh_w <- read_tsv(file_froh_w, show_col_types = FALSE)
df_froh_e <- read_tsv(file_froh_e, show_col_types = FALSE)
df_sift   <- read_tsv(file_sift, show_col_types = FALSE) |>
  filter(Region_type == "all_masked")

# Helper function to merge SIFT and FROH data for both Occidental & Oriental
build_merged <- function(state) {
  df_sift_state <- df_sift |> filter(Allelic_state == state)

  df_w <- inner_join(
    df_froh_w,
    df_sift_state,
    by = c("Sample" = "Sample_ID", "Population" = "Population")
  ) |> mutate(Type = type_mapping[Population], Lineage = "Occidental")

  df_e <- inner_join(
    df_froh_e,
    df_sift_state,
    by = c("Sample" = "Sample_ID", "Population" = "Population")
  ) |> mutate(Type = type_mapping[Population], Lineage = "Oriental")

  bind_rows(df_w, df_e) |>
    mutate(Region = factor(Lineage, levels = c("Occidental", "Oriental")))
}

# ---------------------------------------------------------------------------
# Generic scatter of a load metric vs F_ROH with Region faceting
#   model = "lm"      -> Gaussian linear model (correct for continuous ratios)
#   model = "poisson" -> Poisson GLM (correct only for raw integer counts)
# ---------------------------------------------------------------------------
#   cor_label = NULL  -> ggpubr::stat_cor (raw Pearson p)
#   cor_label = df    -> precomputed Holm-adjusted label (see section 3-8 helpers)
plot_load_vs_froh <- function(df, ycol, ylab, model = "lm", facet = TRUE, cor_label = NULL) {
  fml <- as.formula(paste(ycol, "~ F_ROH"))

  if (model == "poisson") {
    fit <- glm(fml, family = poisson(link = "log"), data = df)
    smooth_layer <- geom_smooth(method = "glm", method.args = list(family = poisson(link = "log")),
                                 color = "gray40", fill = "gray80", alpha = 0.3, linewidth = 1)
  } else {
    fit <- lm(fml, data = df)
    smooth_layer <- geom_smooth(method = "lm",
                                 color = "gray40", fill = "gray80", alpha = 0.3, linewidth = 1)
  }
  print(summary(fit))

  if (is.null(cor_label)) {
    cor_layer <- stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "top",
                          size = 7 / .pt)
  } else {
    cor_layer <- geom_text(data = cor_label, aes(x = mean(range(df$F_ROH)), y = Inf, label = label),
                           parse = TRUE, hjust = 0, vjust = 1.5, size = 7 / .pt, inherit.aes = FALSE)
  }

  p <- ggplot(df, aes(x = F_ROH, y = .data[[ycol]])) +
    geom_point(aes(color = Population, shape = Type),
               size = 2, alpha = 0.8) +
    smooth_layer +
    cor_layer +
    scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE) +
    scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE) +
    theme_classic(base_size = 7) +
    labs(
      x = expression(F[ROH] ~ "(Genome-wide runs of homozygosity proportion)"),
      y = ylab, color = "Population", shape = "Type"
    ) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  if (facet && "Region" %in% names(df)) {
    p <- p + facet_wrap(~ Region, nrow = 1) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )
  }
  return(p)
}

# ---------------------------------------------------------------------------
# 0. Allelic-state configuration
#   All = Ho + He (total load), He = potential load, Ho = realized load.
#   Each state gets its own file-name suffix (_All / _He / _Ho).
# ---------------------------------------------------------------------------
state_config <- list(
  "All" = list(
    ylab_ratio = "Total deleterious burden\n(per synonymous site)",
    ylab_load  = "Total genetic load",
    ylab_count = "Total deleterious mutation count",
    file_tag   = "All",
    ratio_file = "Fig_S_FROH_vs_TotalLoadRatio_perSynSite",
    suffix     = "_All"
  ),
  "He" = list(
    ylab_ratio = "Potential deleterious burden\n(per synonymous site)",
    ylab_load  = "Potential genetic load",
    ylab_count = "Heterozygous deleterious mutation count",
    file_tag   = "Heterozygous",
    ratio_file = "Fig_S_FROH_vs_PotentialLoadRatio_perSynSite",
    suffix     = "_He"
  ),
  "Ho" = list(
    ylab_ratio = "Realized deleterious burden\n(per synonymous site)",
    ylab_load  = "Realized genetic load",
    ylab_count = "Homozygous deleterious mutation count",
    file_tag   = "Homozygous",
    ratio_file = "Fig_S_FROH_vs_RealizedLoadRatio_perSynSite",
    suffix     = "_Ho"
  )
)


# ---------------------------------------------------------------------------
# 2. ROH Boxplot: Cultivar vs Wild with Wilcoxon Test (West East merged)
# ---------------------------------------------------------------------------
df_realized <- build_merged("Ho")

# Subset for West (Occidental) and East (Oriental)
df_realized_w <- df_realized |> filter(Region == "Occidental")
df_realized_e <- df_realized |> filter(Region == "Oriental")

p1_box_w <- ggplot(df_realized_w, aes(x = Type, y = F_ROH)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(shape = Type, color = Population)) +
  scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE) +
  scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE) +
  theme_bw(base_size = 7) +
  labs(title = "Occidental", y = expression(Global ~ F[ROH]), x = "") +
  stat_compare_means(comparisons = list(c("Cultivar", "Wild")), method = "wilcox.test", size = 7 / .pt) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "none"
  )

p1_box_e <- ggplot(df_realized_e, aes(x = Type, y = F_ROH)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(shape = Type, color = Population)) +
  scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE) +
  scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE) +
  theme_bw(base_size = 7) +
  labs(title = "Oriental", y = expression(Global ~ F[ROH]), x = "") +
  stat_compare_means(comparisons = list(c("Cultivar", "Wild")), method = "wilcox.test", size = 7 / .pt) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "none"
  )

box_legend_source <- ggplot(df_realized, aes(x = Type, y = F_ROH)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(shape = Type, color = Population)) +
  scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE) +
  scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_text(face = "bold"), legend.position = "bottom")
box_shared_legend <- ggpubr::get_legend(box_legend_source)

p_box_merge <- (p1_box_w + p1_box_e + plot_layout(ncol = 2, widths = c(1, 1))) /
  wrap_elements(full = box_shared_legend) +
  plot_layout(heights = c(1, 0.3))

save_myplot(p_box_merge, "Fig_S_FROH_boxplot_per_type")

unify_style <- list(
  scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE),
  scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE),
  theme_classic(base_size = 7),
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "none"
  )
)

legend_source <- ggplot(df_realized, aes(x = F_ROH, y = Deleterious_per_synonymous_SNP,
                                          color = Population, shape = Type)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = pop_color, limits = names(pop_color), drop = FALSE) +
  scale_shape_manual(values = type_shape, limits = names(type_shape), drop = FALSE) +
  theme_classic(base_size = 7) +
  theme(legend.title = element_text(face = "bold"), legend.position = "bottom")
shared_legend <- ggpubr::get_legend(legend_source)

df_all_states <- bind_rows(lapply(names(state_config), function(st) {
  build_merged(st) |> mutate(Allelic_state = st)
}))

cor_one <- function(x) {
  t <- cor.test(x$F_ROH, x$Deleterious_per_synonymous_SNP, method = "pearson")
  data.frame(n = nrow(x), R = unname(t$estimate), p = t$p.value)
}

cor_overall <- df_all_states |>
  group_by(Allelic_state, Region) |>
  group_modify(~ cor_one(.x)) |>
  ungroup() |>
  mutate(p_holm = p.adjust(p, method = "holm"))

cor_pop <- df_all_states |>
  group_by(Allelic_state, Region, Population) |>
  filter(n() >= 3) |>
  group_modify(~ cor_one(.x)) |>
  ungroup() |>
  group_by(Allelic_state) |>
  mutate(p_holm = p.adjust(p, method = "holm")) |>
  ungroup()

# plotmath label in stat_cor's style: R = 0.55, p_adj = 1.4e-06
fmt_cor_label <- function(R, p_adj) {
  p_str <- ifelse(p_adj < 2.2e-16, '"<"~"2.2e-16"',
                  paste0('"="~"', vapply(signif(p_adj, 2), format, character(1)), '"'))
  paste0('italic(R)~"="~', signif(R, 2), '*","~italic(p)[adj]~', p_str)
}

cor_overall <- cor_overall |> mutate(label = fmt_cor_label(R, p_holm))
cor_pop     <- cor_pop |> mutate(label = fmt_cor_label(R, p_holm))

print(as.data.frame(cor_overall |> select(-label)), digits = 3)
print(as.data.frame(cor_pop |> select(-label)), digits = 3)
write_tsv(bind_rows(cor_overall |> mutate(Population = "Overall"), cor_pop) |>
            select(Allelic_state, Lineage = Region, Population, n, R, p, p_holm),
          "Fig_ROH_ratio_cor_holm.tsv")

# Per-population labels stacked from the top-right corner, one line per population
pop_label_df <- function(st, region) {
  cor_pop |>
    filter(Allelic_state == st, Region == region) |>
    arrange(Population) |>
    mutate(vjust = 1.5 + 1.3 * (row_number() - 1))
}

plot_state_ratio <- function(state) {
  cfg <- state_config[[state]]
  cat("\n========== Allelic state:", state, "==========\n")

  df_state   <- build_merged(state)
  df_state_w <- df_state |> filter(Region == "Occidental")
  df_state_e <- df_state |> filter(Region == "Oriental")

  # 3. Load Ratio: Occidental & Oriental as side-by-side panels
  p_ratio_w <- plot_load_vs_froh(
    df_state_w,
    ycol = "Deleterious_per_synonymous_SNP",
    ylab = cfg$ylab_ratio,
    model = "lm",
    facet = FALSE,
    cor_label = cor_overall |> filter(Allelic_state == state, Region == "Occidental")
  ) + labs(title = "Occidental")

  p_ratio_e <- plot_load_vs_froh(
    df_state_e,
    ycol = "Deleterious_per_synonymous_SNP",
    ylab = cfg$ylab_ratio,
    model = "lm",
    facet = FALSE,
    cor_label = cor_overall |> filter(Allelic_state == state, Region == "Oriental")
  ) + labs(title = "Oriental")

  # Grouped Regression by Population (panels kept; standalone figure disabled)
  p3_scatter_w <- ggplot(df_state_w, aes(x = F_ROH, y = Deleterious_per_synonymous_SNP, color = Population)) +
    geom_point(aes(shape = Type), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_shape_manual(values = type_shape, limits = c("Cultivar", "Rootstock", "Wild"), drop = FALSE) +
    scale_color_manual(values = pop_color) +
    theme_bw(base_size = 7) +
    labs(
      title = "Occidental",
      x = expression(F[ROH] ~ "(Genome-wide runs of homozygosity proportion)"),
      y = cfg$ylab_ratio
    ) +
    geom_text(data = pop_label_df(state, "Occidental"), aes(x = Inf, y = Inf, label = label, color = Population, vjust = vjust),
              parse = TRUE, hjust = 1.02, size = 7 / .pt, show.legend = FALSE, inherit.aes = FALSE) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  p3_scatter_e <- ggplot(df_state_e, aes(x = F_ROH, y = Deleterious_per_synonymous_SNP, color = Population)) +
    geom_point(aes(shape = Type), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_shape_manual(values = type_shape, limits = c("Cultivar", "Rootstock", "Wild"), drop = FALSE) +
    scale_color_manual(values = pop_color) +
    theme_bw(base_size = 7) +
    labs(
      title = "Oriental",
      x = expression(F[ROH] ~ "(Genome-wide runs of homozygosity proportion)"),
      y = cfg$ylab_ratio
    ) +
    geom_text(data = pop_label_df(state, "Oriental"), aes(x = Inf, y = Inf, label = label, color = Population, vjust = vjust),
              parse = TRUE, hjust = 1.02, size = 7 / .pt, show.legend = FALSE, inherit.aes = FALSE) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

# Combined Plot (Overall Ratio + Per-Population Ratio)
  p_ratio_w_u    <- p_ratio_w + unify_style
  p_ratio_e_u    <- p_ratio_e + unify_style
  p3_scatter_w_u <- p3_scatter_w + unify_style
  p3_scatter_e_u <- p3_scatter_e + unify_style

  panels_grid <- (p_ratio_w_u + p_ratio_e_u) / (p3_scatter_w_u + p3_scatter_e_u)
  ratio_combined_plot <- panels_grid / wrap_elements(full = shared_legend) +
    plot_layout(heights = c(1, 1, 0.35))

  save_myplot(ratio_combined_plot, paste0("Fig_ROH_ratio_overall+per_population", cfg$suffix), w = 180, h = 150)

  # Unified overall + per-population panels, reused by the All/He/Ho overview (section 8)
  return(list(
    ratio_w = p_ratio_w_u, ratio_e = p_ratio_e_u,
    pop_w = p3_scatter_w_u, pop_e = p3_scatter_e_u
  ))
}

ratio_panels <- lapply(names(state_config), plot_state_ratio)
names(ratio_panels) <- names(state_config)

state_rows <- unlist(lapply(ratio_panels, function(x) {
  list(x$ratio_w + x$ratio_e, x$pop_w + x$pop_e)
}), recursive = FALSE)
state_grid <- wrap_plots(state_rows, ncol = 1)
state_combined_plot <- state_grid / wrap_elements(full = shared_legend) +
  plot_layout(heights = c(rep(1, length(state_rows)), 0.35))

save_myplot(state_combined_plot, "Fig_ROH_ratio_All+He+Ho", w = 180, h = 400)

cat("DONE!\n")
