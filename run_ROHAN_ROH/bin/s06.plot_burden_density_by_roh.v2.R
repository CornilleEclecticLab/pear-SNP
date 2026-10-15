library(tidyverse)
library(ggpubr)
library(patchwork)

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# --- config -
in_dir      <- "../output/s05.genetic_burden_by_roh"
ref_genomes <- c("PPY", "PCOM")
color_file  <- "../input/population_color.tsv"   # optional; falls back if absent

denom_col <- "Genotyped_synonymous_sites"

save_myplot <- function(p, prefix, w = 180, h = 120, dpi = 600) {
  wi <- w / 25.4; hi <- h / 25.4
  ggsave(paste0(prefix, ".pdf"), p, width = wi, height = hi, units = "in")
  ggsave(paste0(prefix, ".png"), p, width = wi, height = hi, units = "in", dpi = dpi)
}

type_mapping <- c(
  "ussu"="Wild","pash"="Wild","cauc"="Wild","pyra"="Wild",
  "betu"="Rootstock",
  "comm_Dessert"="Cultivar","comm_Perry"="Cultivar","pyri_JP"="Cultivar",
  "Sand_CN"="Cultivar","Sand_CN-SE"="Cultivar","Sand_CN-SW"="Cultivar","White"="Cultivar"
)
type_shape <- c("Wild"=17, "Rootstock"=7, "Cultivar"=4)
ref_label  <- c("PPY"="Oriental", "PCOM"="Occidental")

pt_size <- function(pt) pt / (72.27 / 25.4)

# --- read numerator, denominator, tag Ref -----------------------------------
# The burden file's `Type` column is zygosity; rename so it does not clash with
# the population `Type` (Wild/Rootstock/Cultivar).
read_ref <- function(ref) {
  bfile <- file.path(in_dir, paste0(ref, "_mutation_burden_by_ROH.txt"))
  dfile <- file.path(in_dir, paste0(ref, "_ROH_denominators.txt"))
  if (!file.exists(bfile) || !file.exists(dfile)) {
    message("Skipping ", ref, " (missing output files)"); return(NULL)
  }
  burden <- read_tsv(bfile, show_col_types = FALSE) |> rename(Zygosity = Type)
  
  denom  <- read_tsv(dfile, show_col_types = FALSE)
  left_join(burden, denom, by = c("Individual", "Population", "Region")) |>
    mutate(Ref = ref)
}

df <- map_dfr(ref_genomes, read_ref)

stopifnot(nrow(df) > 0)

# --- exclude individuals with zero ROH --------------------------------
# Drop both region rows so the paired ROH vs nonROH comparison stays balanced.
roh_sites <- df |>
  filter(Region == "ROH") |>
  group_by(Ref, Individual) |>
  summarise(roh_denom = sum(.data[[denom_col]], na.rm = TRUE), .groups = "drop")

drop_ind <- filter(roh_sites, roh_denom == 0)
if (nrow(drop_ind) > 0) {
  message("Excluding ", nrow(drop_ind), " individual(s) with zero ROH: ",
          paste(unique(drop_ind$Individual), collapse = ", "))
}

df <- semi_join(df, filter(roh_sites, roh_denom > 0),
                by = c("Ref", "Individual"))
stopifnot(nrow(df) > 0)

# --- density -----------
df <- df |>
  mutate(
    denom    = .data[[denom_col]],
    density  = Deleterious_allele_count / denom,
    Region   = factor(Region, levels = c("nonROH", "ROH")),
    # input labels -> display labels (input files unchanged)
    Zygosity = factor(Zygosity, levels = c("Allele", "Homozygous", "Heterozygous"),
                      labels = c("Total", "Homozygous(realized)", "Heterozygous(potential)")),
    Type     = type_mapping[Population]
  ) |>
  filter(is.finite(density))

# --- per-individual wide table + correct paired test -------------------------
wide <- df |>
  select(Ref, Individual, Population, Type, Zygosity, Region, density) |>
  pivot_wider(names_from = Region, values_from = density) |>
  mutate(fold_ROH_over_nonROH = ROH / nonROH)

paired_stats <- wide |>
  group_by(Ref, Zygosity) |>
  summarise(
    n_pairs       = sum(is.finite(ROH) & is.finite(nonROH)),
    median_ROH    = median(ROH, na.rm = TRUE),
    median_nonROH = median(nonROH, na.rm = TRUE),
    p_wilcox      = if (sum(is.finite(ROH) & is.finite(nonROH)) > 1)
                      wilcox.test(ROH, nonROH, paired = TRUE)$p.value else NA_real_,
    .groups = "drop"
  ) |>
  mutate(median_fold = median_ROH / median_nonROH,
         # Holm tests shown in the figure
         p_holm = p.adjust(p_wilcox, method = "holm"))

write_tsv(wide,         "s06_density_by_ROH_per_individual.tsv")
write_tsv(paired_stats, "s06_density_by_ROH_paired_stats.tsv")
print(paired_stats)

# --- plot: paired ROH vs nonROH
pop_color <- NULL
if (file.exists(color_file)) {
  cc <- read_tsv(color_file, show_col_types = FALSE) |> distinct(Population, Color)
  pop_color <- setNames(cc$Color, cc$Population)
}

fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 1e-4, sprintf("p.adj=%.1e", p), sprintf("p.adj=%.4f", p)))

make_plot <- function(ref) {
  dat <- filter(df, Ref == ref)
  if (nrow(dat) == 0) return(NULL)

  lab <- dat |>
    group_by(Zygosity) |>
    summarise(y = max(density, na.rm = TRUE) * 1.05, .groups = "drop") |>
    left_join(filter(paired_stats, Ref == ref), by = "Zygosity") |>
    mutate(label = fmt_p(p_holm), x = 1.5)

  p <- ggplot(dat, aes(x = Region, y = density)) +
    geom_line(aes(group = Individual, color = Population),
              alpha = 0.35, linewidth = 0.4) +
    geom_boxplot(aes(group = Region), width = 0.35, outlier.shape = NA,
                 fill = NA, color = "gray30") +
    geom_point(aes(color = Population, shape = Type), size = 2.5,
               alpha = 0.85, stroke = 1) +
    geom_text(data = lab, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = pt_size(7)) +
    facet_wrap(~ Zygosity, scales = "free_y") +
    scale_shape_manual(values = type_shape) +
    theme_classic(base_size = 7) +
    labs(x = NULL, y = "Deleterious burden\n(per synonymous site)",
         title = paste0(ref_label[ref], ": Deleterious burden, ROH vs non-ROH"),
         subtitle = "Paired Wilcoxon within individual, Holm-adjusted",
         color = "Population", shape = "Type") +
    theme(axis.text = element_text(color = "black"),
          legend.position = "right",
          legend.title = element_text(face = "bold"))
  if (!is.null(pop_color)) p <- p + scale_color_manual(values = pop_color)
  p
}


plot_order <- c("PCOM", "PPY")
plots <- map(plot_order, make_plot) |> set_names(plot_order) |> compact()
stopifnot(length(plots) > 0)

combined <- wrap_plots(plots, ncol = 1)
save_myplot(combined, "s06.Fig_S_load_density_ROH_vs_nonROH", w = 180, h = 180)

cat("DONE!\n")