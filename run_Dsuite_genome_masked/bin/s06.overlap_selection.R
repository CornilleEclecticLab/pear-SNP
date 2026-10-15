library(GenomicRanges)
library(ggrepel)
library(tidyverse)

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


# COnfig paths ----
OUT_DIR <- "../output/Dinvestigate"
GFF_FILE_E <- "../../run_OmegaPlus/input/GWHBAOS00000000.chr_list.gene.slop1kb.gff"
GFF_FILE_W <- "../../run_OmegaPlus_refer_on_communis/input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop1kb.gff"
SELECTION_GENE_FILES_E <- "../../run_OmegaPlus/output/s08.summary_table/any_by_pop.positive_selection_genes_method.*.txt"
SELECTION_GENE_FILES_W <- "../../run_OmegaPlus_refer_on_communis/output/s08.summary_table/any_by_pop.positive_selection_genes_method.*.txt"
MAP_FILE <- "../input/chromosome_accession_ID_map.txt"

PLOT_DIR <- "../plots"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_DIR <- "../output/Dinvestigate_summary"
dir.create(SUMMARY_DIR, showWarnings = FALSE, recursive = TRUE)

# Recognize the pop ID from filename
known_pop_ID <- c("comm_Dessert", "comm_Perry", "cauc", "pyra", "betu", "pash",
                   # longest-first so "Sand_CN-SE" matches before "Sand_CN"
                  "pyri_JP", "Sand_CN-SE", "Sand_CN-SW", "Sand_CN", "ussu", "White") 
pattern_dict <- paste(known_pop_ID[order(-nchar(known_pop_ID))], collapse = "|")

# Identify trios
all_files    <- list.files(OUT_DIR, pattern = "localFstats__50_25\\.txt$",
                           full.names = TRUE, recursive = TRUE)
trio_names   <- str_match(basename(all_files), "^[^.]+\\.(.*)_localFstats__50_25\\.txt$")[, 2]
unique_trios <- unique(trio_names[!is.na(trio_names)])
parsed_trios <- setNames(
  map_chr(str_extract_all(unique_trios, pattern_dict), ~ paste(.x, collapse = ", ")),
  unique_trios
)

# Run for Occidental or Oriental lineage?
if (any(grepl("comm_Dessert|comm_Perry|cauc|pyra", unique_trios))) {
  GFF_FILE <- GFF_FILE_W 
  SELECTION_GENE_FILES <- SELECTION_GENE_FILES_W
  # MASK_PASS_BED <- MASK_PASS_BED_W
  LINEAGE <- "west"
} else {
  GFF_FILE <- GFF_FILE_E 
  SELECTION_GENE_FILES <- SELECTION_GENE_FILES_E
  LINEAGE <- "east"
}


# Statistic-specific outlier definitions
FD_TOP_Q   <- 0.95       # f_d: one-tailed top 5%
FDM_TOP_Q  <- 0.975      # f_dM: two-tailed  bottom 2.5% & top 2.5%

# Load GFF and selection genes
chr_map <- if (file.exists(MAP_FILE)) {
  read_tsv(MAP_FILE, col_names = c("chr_accession", "Chr_ID"), show_col_types = FALSE)
} else NULL

# Using chr ID for ploting rather chr accession
add_chr_label <- function(df, chr_map) {
  if (is.null(chr_map)) {
    df$chr_label <- df$chr  # west: label == native Chr1, Chr2 ...
  } else {
    lut <- setNames(chr_map$Chr_ID, chr_map$chr_accession)
    df$chr_label <- coalesce(unname(lut[df$chr]), df$chr)  # back to native if unmapped
  }
  df
}

# GR for a set of gene IDs
get_gr_from_ids <- function(ids, gff_df) {
  df_sub <- gff_df %>% filter(gene_id %in% ids)
  if (nrow(df_sub) == 0) return(NULL)
  GRanges(df_sub$chr, IRanges(df_sub$start, df_sub$end), gene_id = df_sub$gene_id)
}

# GR of genes selected in p1p2p3 trio of `pops` (intersection).  
# highlight  P2 n P3 and P1 n P3.
# if input P1, P2, P3 as pops, return GR of genes selected in P1nP2 or P1nP3 or P2nP3, or P1nP2nP3
intersect_selection_gr <- function(pops, selection_df, gff_df) {
  pops <- intersect(pops, unique(selection_df$pop))
  if (length(pops) < 2) return(NULL)
  ids <- selection_df %>%
    filter(pop %in% pops) %>%
    distinct(gene_id, pop) %>%
    count(gene_id) %>%
    filter(n >= 2) %>%
    pull(gene_id)
  cat(sprintf("Shared %d genes in {%s}\n", length(ids), paste(pops, collapse = ",")))
  get_gr_from_ids(ids, gff_df)
}

# Annotate each window with overlap hit + comma-joined gene IDs for a GRanges.
# Uses df$chr (NATIVE, Chr1 in west, ChrAccession in east).
# NULL gr -> all FALSE.
annotate_overlap <- function(df, gr) {
  hit   <- rep(FALSE, nrow(df))
  genes <- rep(NA_character_, nrow(df))
  if (!is.null(gr) && nrow(df) > 0) {
    win_gr <- GRanges(df$chr, IRanges(df$windowStart, df$windowEnd))
    ov <- findOverlaps(win_gr, gr)
    if (length(ov) > 0) {
      hit[queryHits(ov)] <- TRUE
      g <- as.data.frame(ov) %>%
        mutate(gene_id = gr$gene_id[subjectHits]) %>%
        group_by(queryHits) %>%
        summarise(s = paste(unique(gene_id), collapse = ","), .groups = "drop")
      genes[g$queryHits] <- g$s
    }
  }
  list(hit = hit, genes = genes)
}

# Load GFF (native chr coords)
prefix <- if (LINEAGE == "east") {
  "Accession="
} else if (LINEAGE == "west") {
  "Name="
}

gff_df <- read_tsv(
  GFF_FILE, 
  col_names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
  comment = "#", 
  show_col_types = FALSE
) %>%
  filter(type == "gene") %>%
  mutate(gene_id = str_match(attributes, paste0(prefix, "([^;]+)"))[, 2]) %>%
  filter(!is.na(gene_id))


# Load positive selection candidate genes
selection_df <- Sys.glob(SELECTION_GENE_FILES) %>%
  map_dfr(read_tsv, col_names = c("gene_id", "pop", "method"), show_col_types = FALSE)

selection_gene_ids <- unique(selection_df$gene_id)
selection_gff <- gff_df %>% filter(gene_id %in% selection_gene_ids)
selection_gr  <- GRanges(selection_gff$chr, IRanges(selection_gff$start, selection_gff$end),
                         gene_id = selection_gff$gene_id)

add_cumulative_pos <- function(df, gap = 2) {
  chr_off <- df %>%
    group_by(chr_num) %>%
    summarise(chr_max = max(mid_pos_Mb, na.rm = TRUE), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(offset = lag(cumsum(chr_max + gap), default = 0))
  df %>%
    left_join(chr_off, by = "chr_num") %>%
    mutate(x_cum = mid_pos_Mb + offset)
}


# Per-population selection GR
pop_selection_gr <- selection_df %>%
  inner_join(gff_df %>% select(gene_id, chr, start, end), by = "gene_id") %>%
  group_split(pop) %>%
  setNames(map_chr(., ~ .x$pop[1])) %>%
  map(~ GRanges(.x$chr, IRanges(.x$start, .x$end), gene_id = .x$gene_id))





# ---- overlap f_d ----------------------------------------------------------
# f_d is one-tailed (positive = gene flow P2<->P3)
# highlight genes selected in both P2 and P3.
walk(unique_trios, function(trio) {
  color_map     <- c(east = "#FFD3D3", west = "#DDD6FF")
  lineage_color <- ifelse(str_detect(trio, "comm|cauc|pyra"), color_map["west"], color_map["east"])
  trio_files    <- all_files[str_detect(basename(all_files), paste0("\\.", trio, "_localFstats"))]

  df_fd <- trio_files %>% map_dfr(~ read_tsv(.x, show_col_types = FALSE)) %>%
    add_chr_label(chr_map) %>%                       # adds chr_label; chr stays native
    mutate(mid_pos_Mb = (windowStart + windowEnd) / 2 / 1e6,
           chr_num    = as.numeric(str_extract(chr_label, "\\d+"))) %>%
    filter(!is.na(chr_num), !is.na(f_d)) %>%
    mutate(chr_label = fct_reorder(chr_label, chr_num)) %>%
    arrange(chr_num, mid_pos_Mb) %>%
    add_cumulative_pos()

  # If too many negative f_d values, transfer to 0, otherwise remove them, same for f_d values > 1
  if (sum(df_fd$f_d < -1, na.rm=TRUE) > 0.01*sum(df_fd$f_d >= 0, na.rm=TRUE)) {
      warning(sprintf("[%s] too many f_d < -1", trio))
  }
  if (sum(df_fd$f_d > 1, na.rm=TRUE) > 0.01*sum(df_fd$f_d >= 0, na.rm=TRUE)) {
      warning(sprintf("[%s] too many f_d > 1", trio))
  }
  
  df_fd <- df_fd %>%
      mutate(f_d = ifelse(f_d < 0, 0, f_d)) %>%
      filter(f_d <= 1)

  # P2 n P3 selection genes
  trio_pops <- str_extract_all(trio, pattern_dict)[[1]]      # P1, P2, P3
  p2p3_gr   <- if (length(trio_pops) >= 3) 
    intersect_selection_gr(c(trio_pops[2], trio_pops[3]), selection_df, gff_df) else NULL

  ann <- annotate_overlap(df_fd, p2p3_gr)
  df_fd <- df_fd %>% mutate(
    is_outlier        = f_d >= quantile(f_d, FD_TOP_Q),
    overlaps_target   = ann$hit,
    overlapping_genes = ann$genes
  )

  # # f_d plot
  # axis breaks at chromosome centers
  axis_df <- df_fd %>%
    group_by(chr_label, chr_num) %>%
    summarise(center = (min(x_cum) + max(x_cum)) / 2, .groups = "drop") %>%
    arrange(chr_num)

  color_palette <- c("0" = "grey", "1" = unname(lineage_color))
  p <- ggplot(df_fd, aes(x_cum, f_d)) +
        geom_point(aes(color = as.factor(chr_num %% 2)), alpha = 0.4, size = 0.8, stroke = 0) +
        geom_hline(yintercept = quantile(df_fd$f_d, FD_TOP_Q), color = "black", linetype = "dashed", linewidth = 0.5) +
        scale_color_manual(values = color_palette, guide = "none")

  df_overlap <- df_fd %>% filter(overlaps_target)
  if (nrow(df_overlap) > 0) {
    p <- p + geom_point(data = df_overlap, aes(x_cum, f_d), color = "blue", alpha = 0.8, size = 1.5, shape = 17)
  }

  df_both <- df_fd %>% filter(is_outlier & overlaps_target)
  if (nrow(df_both) > 0) {
    df_labels <- df_both %>%
      separate_rows(overlapping_genes, sep = ",") %>%
      group_by(overlapping_genes) %>%
      slice_max(f_d, n = 1, with_ties = FALSE) %>%
      ungroup()

    p <- p +
      geom_point(data = df_both, aes(x_cum, f_d), color = "#E63323", alpha = 1, size = 1.8, shape = 18) +
      geom_text_repel(
        data = df_labels, aes(x_cum, f_d, label = overlapping_genes),
        size = 2.5, color = "#E63323",
        box.padding = 0.4, point.padding = 0.2,
        min.segment.length = 0, segment.size = 0.2, segment.color = "grey60",
        direction = "both", max.overlaps = Inf, seed = 1)
  }

  p <- p +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_label,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 7) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.margin = margin(t = 14, r = 12, b = 6, l = 12)) +
    labs(title = paste("f_d on Trios:", parsed_trios[[trio]], "(P1,P2,P3)"),
         subtitle = "Positive = gene flow between P2 P3. Blue triangles/red diamonds overlap genes selected in both P2 and P3.",
         x = "Chromosome", y = expression(italic(f)[d]))

  save_myplot(p, file.path(PLOT_DIR, paste0("Manhattan_fd_", trio)), w = 420, h = 120)
})

# ---- overlap f_dM ---------------------------------------------------------
# f_dM is two-tailed: positive = P2<->P3, negative = P1<->P3. 
# Red diamonds are sign-matched: 
#   positive outlier overlaps with P2 n P3
#   negative outlier overlaps with P1 n P3.
walk(unique_trios, function(trio) {
  color_map     <- c(east = "#FFD3D3", west = "#DDD6FF")
  lineage_color <- ifelse(str_detect(trio, "comm|cauc|pyra"), color_map["west"], color_map["east"])
  trio_files    <- all_files[str_detect(basename(all_files), paste0("\\.", trio, "_localFstats"))]

  df <- trio_files %>% map_dfr(~ read_tsv(.x, show_col_types = FALSE)) %>%
    add_chr_label(chr_map) %>%
    mutate(mid_pos_Mb = (windowStart + windowEnd) / 2 / 1e6,
           chr_num    = as.numeric(str_extract(chr_label, "\\d+"))) %>%
    filter(!is.na(chr_num), !is.na(f_dM)) %>%
    mutate(chr_label = fct_reorder(chr_label, chr_num)) %>%
    arrange(chr_num, mid_pos_Mb) %>%
    add_cumulative_pos()
    
  # f_dM bounded [-1,1]
  # If too many f_dM <-1 values, transfer to 0, otherwise remove them, same for f_dM values > 1
  if (sum(df$f_dM < -1, na.rm=TRUE) > 0.01*sum(df$f_dM >= 0, na.rm=TRUE)) {
      warning(sprintf("[%s] too many f_dM < -1", trio))
  }
  if (sum(df$f_dM > 1, na.rm=TRUE) > 0.01*sum(df$f_dM >= 0, na.rm=TRUE)) {
      warning(sprintf("[%s] too many f_dM > 1", trio))
  }

  df <- df %>%
      filter(f_dM >= -1) %>%
      filter(f_dM <= 1)


  axis_df <- df %>%
    group_by(chr_label, chr_num) %>%
    summarise(center = (min(x_cum) + max(x_cum)) / 2, .groups = "drop") %>%
    arrange(chr_num)
  
  q_low  <- quantile(df$f_dM, 1 - FDM_TOP_Q) 
  q_high <- quantile(df$f_dM, FDM_TOP_Q)

  trio_pops <- str_extract_all(trio, pattern_dict)[[1]]   # P1, P2, P3
  p2p3_gr <- if (length(trio_pops) >= 3)
    intersect_selection_gr(c(trio_pops[2], trio_pops[3]), selection_df, gff_df) else NULL
  p1p3_gr <- if (length(trio_pops) >= 3)
    intersect_selection_gr(c(trio_pops[1], trio_pops[3]), selection_df, gff_df) else NULL

  ann_pos <- annotate_overlap(df, p2p3_gr)   # for positive outliers
  ann_neg <- annotate_overlap(df, p1p3_gr)   # for negative outliers

  df <- df %>% mutate(
    is_pos = f_dM >= q_high,
    is_neg = f_dM <= q_low,
    ov_pos = ann_pos$hit, g_pos = ann_pos$genes,
    ov_neg = ann_neg$hit, g_neg = ann_neg$genes,
    overlaps_target = ov_pos | ov_neg,                       # blue: any relevant set
    is_diamond      = (is_pos & ov_pos) | (is_neg & ov_neg), # red: sign-matched
    overlapping_genes = case_when(
      is_pos & ov_pos ~ g_pos,
      is_neg & ov_neg ~ g_neg,
      TRUE ~ NA_character_
    )
  )

  # color_palette <- c("0" = "grey", "1" = unname(lineage_color))
  # p <- ggplot(df, aes(mid_pos_Mb, f_dM)) +
  #       geom_point(aes(color = as.factor(chr_num %% 2)), alpha = 0.4, size = 0.8, stroke = 0) +
  #       # geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "#E63323", linewidth = 0.6) +
  #       geom_hline(yintercept = q_low, color = "black", linetype = "dashed", linewidth = 0.5) +
  #       geom_hline(yintercept = q_high, color = "black", linetype = "dashed", linewidth = 0.5) +
  #       scale_color_manual(values = color_palette, guide = "none")

  #   df_overlap <- df %>% filter(overlaps_target)
  # if (nrow(df_overlap) > 0)
  #   p <- p + geom_point(data = df_overlap, color = "blue", alpha = 0.8, size = 1.2, shape = 17)

  # df_both <- df %>% filter(is_diamond)
  # if (nrow(df_both) > 0) {
  #   df_labels <- df_both %>%
  #     separate_rows(overlapping_genes, sep = ",") %>%
  #     group_by(overlapping_genes) %>%
  #     slice_max(abs(f_dM), n = 1, with_ties = FALSE) %>%
  #     ungroup()

  #   p <- p +
  #     geom_point(data = df_both, color = "#E63323", alpha = 1, size = 1.8, shape = 18) +
  #     geom_text_repel(data = df_labels, aes(label = overlapping_genes),
  #               size = 2.5, color = "#E63323", nudge_y = 0.05, max.overlaps = Inf)
  # }

  # p <- p +
  #   facet_grid(. ~ chr_label, scales = "free_x", space = "free_x", switch = "x") +
  #   theme_minimal(base_size = 7) +
  #   theme(panel.spacing.x = unit(0, "lines"),
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         strip.background = element_blank(),
  #         strip.placement = "outside",
  #         strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
  #         panel.grid.minor = element_blank(),
  #         panel.grid.major.x = element_blank(),
  #         plot.title = element_text(face = "bold")) +
  #   labs(title = paste("f_dM on Trios:", parsed_trios[[trio]], "(P1,P2,P3)"),
  #        subtitle = "Positive: gene flow between P2 P3, Negative: gene flow between P1 P3.\nRed diamonds: sign-matched outlier overlapping genes selected in both P2&P3 (positive) or P1&P3 (negative).",
  #        x = "Chromosome", y = expression(italic(f)[dM]))
  # save_myplot(p, file.path(PLOT_DIR, paste0("Manhattan_fdM_", trio)), w = 420, h = 120)
  
  color_palette <- c("0" = "grey", "1" = unname(lineage_color))
  p <- ggplot(df, aes(x_cum, f_dM)) +
        geom_point(aes(color = as.factor(chr_num %% 2)), alpha = 0.4, size = 0.8, stroke = 0) +
        geom_hline(yintercept = q_low,  color = "black", linetype = "dashed", linewidth = 0.5) +
        geom_hline(yintercept = q_high, color = "black", linetype = "dashed", linewidth = 0.5) +
        scale_color_manual(values = color_palette, guide = "none")

  df_overlap <- df %>% filter(overlaps_target)
  if (nrow(df_overlap) > 0)
    p <- p + geom_point(data = df_overlap, aes(x_cum, f_dM), color = "blue", alpha = 0.8, size = 1.5, shape = 17)

  df_both <- df %>% filter(is_diamond)
  if (nrow(df_both) > 0) {
    df_labels <- df_both %>%
      separate_rows(overlapping_genes, sep = ",") %>%
      group_by(overlapping_genes) %>%
      slice_max(abs(f_dM), n = 1, with_ties = FALSE) %>%
      ungroup()

    p <- p +
      geom_point(data = df_both, aes(x_cum, f_dM), color = "#E63323", alpha = 1, size = 1.8, shape = 18) +
      geom_text_repel(
        data = df_labels, aes(x_cum, f_dM, label = overlapping_genes),
        size = 2.5, color = "#E63323",
        box.padding = 0.4, point.padding = 0.2,
        min.segment.length = 0, segment.size = 0.2, segment.color = "grey60",
        direction = "both", max.overlaps = Inf, seed = 1)
  }

  p <- p +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_label,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.20, 0.30))) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 7) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.margin = margin(t = 14, r = 12, b = 6, l = 12)) +
    labs(title = paste("f_dM on Trios:", parsed_trios[[trio]], "(P1,P2,P3)"),
         subtitle = "Positive: gene flow between P2 P3, Negative: gene flow between P1 P3.\nRed diamonds: sign-matched outlier overlapping genes selected in both P2&P3 (positive) or P1&P3 (negative).",
         x = "Chromosome", y = expression(italic(f)[dM]))

  save_myplot(p, file.path(PLOT_DIR, paste0("Manhattan_fdM_", trio)), w = 420, h = 120)
})

# Per-trio LocalFstats summary table ------
dir.create(SUMMARY_DIR, showWarnings = FALSE, recursive = TRUE)

walk(unique_trios, function(trio) {
  trio_files <- all_files[str_detect(basename(all_files), paste0("\\.", trio, "_localFstats"))]
  df <- trio_files %>% map_dfr(~ read_tsv(.x, show_col_types = FALSE))
  if (nrow(df) == 0) return(invisible(NULL))

  trio_pops <- str_extract_all(trio, pattern_dict)[[1]]   # P1, P2, P3
  p2p3_gr <- if (length(trio_pops) >= 3)
    intersect_selection_gr(c(trio_pops[2], trio_pops[3]), selection_df, gff_df) else NULL
  p1p3_gr <- if (length(trio_pops) >= 3)
    intersect_selection_gr(c(trio_pops[1], trio_pops[3]), selection_df, gff_df) else NULL

  # window-level overlap
  ann_p2p3 <- annotate_overlap(df, p2p3_gr)   # for positive f_d and f_dM
  ann_p1p3 <- annotate_overlap(df, p1p3_gr)   # for negative f_dM

  # f_d outlier (single tail)
  valid_fd   <- !is.na(df$f_d) & df$f_d <= 1
  fd_clamped <- pmax(df$f_d, 0)
  q_fd       <- if (any(valid_fd)) unname(quantile(fd_clamped[valid_fd], FD_TOP_Q)) else NA_real_
  is_fd_out  <- valid_fd & !is.na(q_fd) & fd_clamped >= q_fd

  # f_dM outlier (two tails)
  valid_fdm  <- !is.na(df$f_dM) & df$f_dM >= -1 & df$f_dM <= 1
  q_lo <- if (any(valid_fdm)) unname(quantile(df$f_dM[valid_fdm], 1 - FDM_TOP_Q)) else NA_real_
  q_hi <- if (any(valid_fdm)) unname(quantile(df$f_dM[valid_fdm], FDM_TOP_Q))     else NA_real_
  is_fdm_top <- valid_fdm & !is.na(q_hi) & df$f_dM >= q_hi
  is_fdm_bot <- valid_fdm & !is.na(q_lo) & df$f_dM <= q_lo

  df_out <- df %>% mutate(
  "{paste0('If_top_', 100 - FD_TOP_Q*100, '_percent_f_d_outlier')}" := is_fd_out,
  f_d_outlier_overlap_selection_gene_ID = ifelse(is_fd_out & ann_p2p3$hit, ann_p2p3$genes, NA_character_),
  
  "{paste0('If_top_', 100 - FDM_TOP_Q*100, '_percent_f_dM_outlier')}" := is_fdm_top,
  f_dM_top_outlier_overlap_selection_gene_ID = ifelse(is_fdm_top & ann_p2p3$hit, ann_p2p3$genes, NA_character_),
  
  "{paste0('If_bottom_', 100 - FDM_TOP_Q*100, '_percent_f_dM_outlier')}" := is_fdm_bot,
  f_dM_bottom_outlier_overlap_selection_gene_ID = ifelse(is_fdm_bot & ann_p1p3$hit, ann_p1p3$genes, NA_character_)
  )

  # Record trio identity so P1/P2/P3 are explicit in the table (trio_pops is c(P1,P2,P3);
  # out-of-range index returns NA, so trios with <3 recognised pops stay NA-safe).
  df_out <- df_out %>%
    mutate(Trio = paste0(trio_pops[1],",",trio_pops[2],",",trio_pops[3]),
           P1 = trio_pops[1], P2 = trio_pops[2], P3 = trio_pops[3],
           .before = 1)

  write_tsv(df_out, file.path(SUMMARY_DIR, paste0(trio, "_LocalFstats_summary.txt")))

  df_outlier <- df_out %>%
    filter(is_fd_out | is_fdm_top | is_fdm_bot)

  write_tsv(df_outlier, file.path(SUMMARY_DIR, paste0(trio, "_LocalFstats_summary.outlier.txt")))
})

cat("Done")