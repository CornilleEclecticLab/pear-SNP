#!/usr/bin/env Rscript

################################################################
# s07: enrichment of introgression outliers in selective sweeps
################################################################

# statistical enrichment of introgression-outlier windows in
# shared se regions, on a mappability-masked universe.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
})

# ---- config
SUMMARY_DIR <- "../output/Dinvestigate_summary"

# sweep BED (positive selection regions) + genmap PASS mask, per lineage
SELECTION_BED_GLOB <- "../input/selection_bed/selection.merge_across_methods.*.merged.bed"
MASK_PASS_BED_E      <- "../input/GWHBAOS00000000.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed"
MASK_PASS_BED_W      <- "../input/PyrusCommunis_BartlettDHv2.0.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed"

FD_TOP_Q      <- 0.95    # f_d : one-tailed top 5%     (match s06)
FDM_TOP_Q     <- 0.975   # f_dM: two-tailed 2.5% each  (match s06)
N_PERM        <- 10000   # min attainable p = 1/(N_PERM+1)
SEED          <- 42
PASS_MIN_FRAC <- 0.8     # window >=80% covered by genmap PASS to enter the universe
N_CORE        <- 4L      # >1 -> mclapply over trios; BH stays global

known_pop_ID <- c("comm_Dessert", "comm_Perry", "cauc", "pyra", "betu", "pash",
                  "pyri_JP", "Sand_CN-SE", "Sand_CN-SW", "Sand_CN", "ussu", "White")
pattern_dict <- paste(known_pop_ID[order(-nchar(known_pop_ID))], collapse = "|")

# ---- discover trios; assign lineage
read_summary  <- function(fp) read_tsv(fp, show_col_types = FALSE, guess_max = Inf)
summary_files <- Sys.glob(file.path(SUMMARY_DIR, "*_LocalFstats_summary.txt"))
stopifnot(length(summary_files) > 0)
unique_trios  <- str_replace(basename(summary_files), "_LocalFstats_summary\\.txt$", "")

if (any(grepl("comm_Dessert|comm_Perry|cauc|pyra", unique_trios))) {
  MASK_PASS_BED <- MASK_PASS_BED_W
  LINEAGE <- "Occidental"
} else {
  MASK_PASS_BED <- MASK_PASS_BED_E
  LINEAGE <- "Oriental"
}

# Pops present in run (Occidentl or Oriental) trios. 
# Pop names are lineage-specific, so filtering
run_pops <- unique(unlist(str_extract_all(unique_trios, pattern_dict)))
stopifnot(length(run_pops) > 0)

# ---- helpers
make_win_gr <- function(df) GRanges(df$chr, IRanges(df$windowStart, df$windowEnd))

pass_frac <- function(win_gr, pass_gr) {  # fraction of each window covered by PASS
  hits <- findOverlaps(win_gr, pass_gr)
  out  <- numeric(length(win_gr))
  if (length(hits) > 0) {
    ov <- width(pintersect(win_gr[queryHits(hits)], pass_gr[subjectHits(hits)]))
    s  <- tapply(ov, queryHits(hits), sum)
    out[as.integer(names(s))] <- s
  }
  out / width(win_gr)
}

load_selection_beds <- function(glob_pattern, keep_pops = NULL) {  # pop is 3rd dot-field
  files <- Sys.glob(glob_pattern)
  stopifnot(length(files) > 0)
  pops <- map_chr(basename(files), ~ str_split(.x, "\\.")[[1]][3])
  if (!is.null(keep_pops)) {                       # keep only this run's (same-lineage) pops
    sel   <- pops %in% keep_pops
    files <- files[sel]; pops <- pops[sel]
  }
  if (length(files) == 0)
    stop("No selection BED matched the run's pops: ", paste(keep_pops, collapse = ", "))
  set_names(map(files, function(f) {
    d <- read_tsv(f, col_names = c("chr", "start", "end"), col_types = "cii",
                  comment = "#", progress = FALSE)
    GenomicRanges::reduce(GRanges(d$chr, IRanges(d$start + 1L, d$end)))   # BED 0- -> 1-based
  }), pops)
}

shared_sweep <- function(pops, gr_list) {
  pops <- intersect(pops, names(gr_list))
  if (length(pops) < 2) return(GRanges())
  Reduce(GenomicRanges::intersect, gr_list[pops])
}

circ_perm <- function(win_gr, x, target_gr, stat = c("count", "mean_diff"),
                      n_perm = N_PERM, seed = SEED) {
  stat <- match.arg(stat)
  empty <- tibble(observed = NA_real_, expected = NA_real_, fold = NA_real_,
                  z = NA_real_, p_perm = NA_real_, n_win_in_sweep = 0L)
  if (length(target_gr) == 0 || length(x) == 0) return(empty)

  ord    <- order(as.character(seqnames(win_gr)), start(win_gr))
  win_gr <- win_gr[ord]; x <- x[ord]
  mask   <- overlapsAny(win_gr, target_gr)
  if (!any(mask) || all(mask)) return(empty)
  idx_by_chr <- split(seq_along(x), as.character(seqnames(win_gr)))

  f <- switch(stat,
    count     = function(v) sum(v[mask]),
    mean_diff = function(v) mean(v[mask]) - mean(v[!mask]))
  obs <- f(x)

  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) {
    v <- x
    for (ix in idx_by_chr) {
      n <- length(ix); if (n < 2) next
      k <- sample.int(n, 1L)
      v[ix] <- x[ix][((seq_len(n) - 1L + k) %% n) + 1L]
    }
    f(v)
  }, numeric(1))

  exp_ <- mean(null); sdn <- sd(null)
  tibble(observed = as.numeric(obs), expected = exp_,
         fold   = if (stat == "count" && exp_ > 0) obs / exp_ else NA_real_,
         z      = if (!is.na(sdn) && sdn > 0) (obs - exp_) / sdn else NA_real_,
         p_perm = (sum(null >= obs) + 1) / (n_perm + 1),
         n_win_in_sweep = sum(mask))
}

thin_windows <- function(win_gr) {
  ord <- order(as.character(seqnames(win_gr)), start(win_gr))
  keep <- logical(length(win_gr)); last_end <- -Inf; last_chr <- ""
  for (i in ord) {
    ch <- as.character(seqnames(win_gr)[i])
    if (ch != last_chr || start(win_gr)[i] > last_end) {
      keep[i] <- TRUE; last_end <- end(win_gr)[i]; last_chr <- ch
    }
  }
  keep
}

# Fisher on LD-thinned windows. 
# Also returns the 2x2 cell counts
fisher_thinned <- function(win_gr, is_out, target_gr) {
  na_row <- tibble(or_fisher = NA_real_, p_fisher = NA_real_, n_win_thinned = NA_integer_,
                   fish_out_in = NA_integer_, fish_out_out = NA_integer_,
                   fish_nonout_in = NA_integer_, fish_nonout_out = NA_integer_)
  if (length(target_gr) == 0 || length(is_out) == 0) return(na_row)

  k        <- thin_windows(win_gr)
  is_out_k <- as.logical(is_out[k]) %in% TRUE          # NA -> FALSE
  hit      <- overlapsAny(win_gr, target_gr)[k]

  a <- sum( is_out_k &  hit)   # outlier,     in sweep
  b <- sum( is_out_k & !hit)   # outlier,     out of sweep
  c <- sum(!is_out_k &  hit)   # non-outlier, in sweep
  d <- sum(!is_out_k & !hit)   # non-outlier, out of sweep
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)   # rows: outlier T/F ; cols: in/out sweep

  fr <- tryCatch(fisher.test(m, alternative = "greater"), error = function(e) NULL)
  tibble(or_fisher      = if (is.null(fr)) NA_real_ else unname(fr$estimate),
         p_fisher       = if (is.null(fr)) NA_real_ else fr$p.value,
         n_win_thinned  = sum(k),
         fish_out_in    = a, fish_out_out    = b,
         fish_nonout_in = c, fish_nonout_out = d)
}

# ---- build mask + sweeps, restrict both to PASS space
pass_gr <- {
  d <- read_tsv(MASK_PASS_BED, col_names = c("chr", "start", "end"),
                col_types = "cii", comment = "#", progress = FALSE)
  GenomicRanges::reduce(GRanges(d$chr, IRanges(d$start + 1L, d$end)))
}
sweep_gr_list <- load_selection_beds(SELECTION_BED_GLOB, keep_pops = run_pops)
# sweep target clipped to PASS, so target and window universe share one space
sweep_gr_list <- map(sweep_gr_list, ~ GenomicRanges::intersect(.x, pass_gr))

chr_dinv <- unique(read_summary(summary_files[1])$chr)
if (mean(seqlevels(pass_gr) %in% chr_dinv) < 0.8)
  stop("MASK_PASS_BED chr names do not match windows.")
walk(names(sweep_gr_list), function(p) {
  if (mean(seqlevels(sweep_gr_list[[p]]) %in% chr_dinv) < 0.8)
    stop("Selection BED chr names for pop '", p, "' do not match windows.")
})

# Outlier-flag columns written by s06)
COL_FD_OUT      <- paste0("If_top_",    100 - FD_TOP_Q  * 100, "_percent_f_d_outlier")
COL_FDM_TOP_OUT <- paste0("If_top_",    100 - FDM_TOP_Q * 100, "_percent_f_dM_outlier")
COL_FDM_BOT_OUT <- paste0("If_bottom_", 100 - FDM_TOP_Q * 100, "_percent_f_dM_outlier")

# ---- trio-level enrichment (three complementary tests per row):
#   count     circ_perm  -> categorical: are outlier windows enriched in sweeps?
#   mean_diff circ_perm  -> continuous : is mean f_d / f_dM higher inside sweeps?
#   fisher_thinned       -> categorical association on LD-thinned windows
enrich_one <- function(trio) {
  fp <- file.path(SUMMARY_DIR, paste0(trio, "_LocalFstats_summary.txt"))
  if (!file.exists(fp)) return(NULL)
  df <- read_summary(fp)
  trio_pops <- str_extract_all(trio, pattern_dict)[[1]]
  if (length(trio_pops) < 3) return(NULL)

  gr_p2p3 <- shared_sweep(c(trio_pops[2], trio_pops[3]), sweep_gr_list)
  gr_p1p3 <- shared_sweep(c(trio_pops[1], trio_pops[3]), sweep_gr_list)

  win_gr  <- make_win_gr(df)
  ok_pass <- pass_frac(win_gr, pass_gr) >= PASS_MIN_FRAC   # universe restricted to PASS
  fd      <- pmax(df$f_d, 0)

  ok_fd  <- ok_pass & !is.na(df$f_d)  & df$f_d  <= 1
  ok_fdm <- ok_pass & !is.na(df$f_dM) & df$f_dM >= -1 & df$f_dM <= 1

  id <- function(metric, tail, target)
    tibble(trio = trio, lineage = LINEAGE,
           P1 = trio_pops[1], P2 = trio_pops[2], P3 = trio_pops[3],
           metric = metric, tail = tail, target = target)

  r1 <- bind_cols(
    id("f_d", "top", "P2nP3"),
    circ_perm(win_gr[ok_fd], df[[COL_FD_OUT]][ok_fd], gr_p2p3, "count",     seed = SEED + 1L),
    circ_perm(win_gr[ok_fd], fd[ok_fd],               gr_p2p3, "mean_diff", seed = SEED + 2L) %>%
      select(cont_obs_diff = observed, cont_exp_diff = expected, cont_z = z, p_perm_cont = p_perm),
    fisher_thinned(win_gr[ok_fd], df[[COL_FD_OUT]][ok_fd], gr_p2p3))

  r2 <- bind_cols(
    id("f_dM", "top", "P2nP3"),
    circ_perm(win_gr[ok_fdm], df[[COL_FDM_TOP_OUT]][ok_fdm], gr_p2p3, "count",     seed = SEED + 3L),
    circ_perm(win_gr[ok_fdm], df$f_dM[ok_fdm],               gr_p2p3, "mean_diff", seed = SEED + 4L) %>%
      select(cont_obs_diff = observed, cont_exp_diff = expected, cont_z = z, p_perm_cont = p_perm),
    fisher_thinned(win_gr[ok_fdm], df[[COL_FDM_TOP_OUT]][ok_fdm], gr_p2p3))

  r3 <- bind_cols(
    id("f_dM", "bottom", "P1nP3"),
    circ_perm(win_gr[ok_fdm], df[[COL_FDM_BOT_OUT]][ok_fdm], gr_p1p3, "count",     seed = SEED + 5L),
    circ_perm(win_gr[ok_fdm], -df$f_dM[ok_fdm],              gr_p1p3, "mean_diff", seed = SEED + 6L) %>%
      select(cont_obs_diff = observed, cont_exp_diff = expected, cont_z = z, p_perm_cont = p_perm),
    fisher_thinned(win_gr[ok_fdm], df[[COL_FDM_BOT_OUT]][ok_fdm], gr_p1p3))

  bind_rows(r1, r2, r3)
}

if (N_CORE > 1L) {
  enrich_raw <- bind_rows(parallel::mclapply(unique_trios, enrich_one, mc.cores = N_CORE))
} else {
  enrich_raw <- map_dfr(unique_trios, enrich_one)
}

# Global BH across every trio x metric x tail
enrich_tbl <- enrich_raw %>%
  mutate(p_perm_BH      = p.adjust(p_perm,      method = "BH"),
         p_perm_cont_BH = p.adjust(p_perm_cont, method = "BH"),
         p_fisher_BH    = p.adjust(p_fisher,    method = "BH")) %>%
  relocate(p_perm_BH,      .after = p_perm) %>%
  relocate(p_perm_cont_BH, .after = p_perm_cont) %>%
  relocate(p_fisher_BH,    .after = p_fisher) %>%
  arrange(p_perm)

write_tsv(enrich_tbl, file.path(SUMMARY_DIR, "sweep_geneflow_enrichment.tsv"))
cat(sprintf("Done: %d rows (%d trios x 3 tests)\nsweep_geneflow_enrichment.tsv\n",
            nrow(enrich_tbl), length(unique_trios)))
