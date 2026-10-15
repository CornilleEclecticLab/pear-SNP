#!/usr/bin/env Rscript

################################################################
# s08.x1 : pi-matched enrichment of introgression outliers in sweeps
################################################################

# Does the sweep ~ introgression-outlier association of s07.x1 survive
# once local nucleotide diversity is controlled for?
#
# Why: Dinvestigate windows hold a fixed number of SNPs (50), so they are
# physically longer where diversity is low. f_d / f_dM are inflated at low
# diversity, and sweeps are low-diversity by definition. The circular
# permutation of s07.x1 keeps sweep clustering but not pi, so here pi (P1, P2,
# P3) enters as an explicit covariate (logistic regression) and as a matching
# variable (pi-matched non-sweep controls).
#
# Outlier flags and targets are the ones of s06 / s07.x1:
#   f_d  top 5%     vs sweeps shared by P2 & P3
#   f_dM top 2.5%   vs sweeps shared by P2 & P3
#   f_dM bottom 2.5% vs sweeps shared by P1 & P3
# Flags are computed per trio on ALL windows (as s06 does), then the PASS
# universe (>= PASS_MIN_FRAC) is applied (as s07.x1 does).
#
# Recombination rate is NOT available for either reference, so it is not in
# the model. Window length (50 SNPs / bp) partly absorbs local SNP density,
# but not recombination; this is logged at run time.
#
# Run one lineage at a time (Occidental and Oriental use different references).
#
# Steps
#   1 export Dinvestigate windows as BED (-> pixy --bed_file)
#   2 pi per window, P1/P2/P3  (--pi_mode aggregate | bed)
#   3 per-window table: covariates, in_sweep, outlier flags
#   4 pi-matched controls + Fisher
#   5 logistic regression; in_sweep tested by U = O - E (score), p from circular permutation of sweeps
#   6 window-length check


suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(tidyverse)
})

# ---- config (defaults; all overridable on the CLI)
FD_TOP_Q        <- 0.95    # f_d : one-tailed top 5%     (match s06)
FDM_TOP_Q       <- 0.975   # f_dM: two-tailed 2.5% each  (match s06)
PASS_MIN_FRAC   <- 0.8     # window >= 80% PASS to enter the universe
SNPS_PER_WINDOW <- 50L     # Dinvestigate window size (localFstats__50_25)
N_PERM          <- 10000   # min attainable p = 1/(N_PERM+1)
SEED            <- 42
N_MATCH         <- 5L
PI_TOL          <- 0.1     # caliper on the natural-log scale
MIN_PIXY_COV    <- 0.8     # aggregate mode: window fraction covered by pixy windows
MAX_PIXY_RATIO  <- 5       # aggregate mode: stop if median pixy win > this x median Dinv win
OUT_DIR_DEFAULT <- "../output/s08.x1.enrichment_pi_matched"

known_pop_ID <- c("comm_Dessert", "comm_Perry", "cauc", "pyra", "betu", "pash",
                  "pyri_JP", "Sand_CN-SE", "Sand_CN-SW", "Sand_CN", "ussu", "White")
pattern_dict <- paste(known_pop_ID[order(-nchar(known_pop_ID))], collapse = "|")
OCC_PATTERN  <- "comm_Dessert|comm_Perry|cauc|pyra"   # same lineage rule as s07.x1

COL_IN  <- "#eb6834"   # in sweep
COL_OUT <- "#2a78d6"   # outside sweep

# ---- CLI
opt_list <- list(
  make_option("--lineage", type = "character",
              help = "Occidental | Oriental (required; checked against trio names)"),
  make_option("--dinv", type = "character",
              help = "Dinvestigate localFstats files: glob(s), comma-separated. Trio parsed from '<chr>.<P1_P2_P3>_localFstats*.txt'"),
  make_option("--sweep_bed", type = "character",
              help = "Sweep BED glob. Several files -> per-pop BEDs, pop = 3rd dot-field (as s07.x1), targets = shared P2&P3 / P1&P3. One file -> used as the target for every test"),
  make_option("--mask_bed", type = "character", help = "genmap PASS mask BED (0-based)"),
  make_option("--genes", type = "character",
              help = "GFF/GFF3 (type == gene used) or gene BED (0-based), for gene density"),
  make_option("--pixy_pi", type = "character",
              help = "pixy *_pi.txt glob(s), comma-separated"),
  make_option("--pop_map", type = "character", default = NULL,
              help = "TSV with header: trio, P1, P2, P3 (pixy pop names). Default: pop names parsed from the trio"),
  make_option("--pi_mode", type = "character", default = "aggregate",
              help = "aggregate (default) | bed  [%default]"),
  make_option("--min_pixy_cov", type = "double", default = MIN_PIXY_COV,
              help = "aggregate mode: min fraction of a window covered by pixy windows [%default]"),
  make_option("--max_pixy_ratio", type = "double", default = MAX_PIXY_RATIO,
              help = "aggregate mode: stop if median pixy window is > this x median Dinvestigate window [%default]"),
  make_option("--allow_coarse_pi", action = "store_true", default = FALSE,
              help = "aggregate mode: only warn (do not stop) when pixy windows are too coarse"),
  make_option("--pass_min_frac", type = "double", default = PASS_MIN_FRAC, help = "[%default]"),
  make_option("--min_sweep_overlap", type = "integer", default = 1L,
              help = "bp of overlap with a sweep to call in_sweep [%default]"),
  make_option("--n_match", type = "integer", default = N_MATCH, help = "[%default]"),
  make_option("--pi_tol", type = "double", default = PI_TOL,
              help = "matching caliper on log pi and log window_length [%default]"),
  make_option("--no_thin", action = "store_true", default = FALSE,
              help = "do NOT LD-thin windows before matching (s07.x1 thins before Fisher)"),
  make_option("--pi_pseudo", type = "double", default = NA,
              help = "added to pi before log; default = half the smallest non-zero pi"),
  make_option("--chrom_sizes", type = "character", default = NULL,
              help = "2-col chr/length (e.g. .fai) for rotation; default = max end seen in inputs"),
  make_option("--n_perm", type = "integer", default = N_PERM, help = "[%default]"),
  make_option("--threads", type = "integer", default = 1L, help = "[%default]"),
  make_option("--seed", type = "integer", default = SEED, help = "[%default]"),
  make_option("--trios", type = "character", default = NULL,
              help = "optional comma-separated subset of trios"),
  make_option("--s06_summary_dir", type = "character", default = NULL,
              help = "optional: cross-check recomputed outlier flags against s06 *_LocalFstats_summary.txt"),
  make_option("--out_dir", type = "character", default = OUT_DIR_DEFAULT, help = "[%default]"),
  make_option("--tag", type = "character", default = NULL,
              help = "run tag for output names; default <lineage>.<pi_mode>"),
  make_option("--export_windows_only", action = "store_true", default = FALSE,
              help = "run step 1 only (window BED for pixy), then exit"),
  make_option("--dry_run", action = "store_true", default = FALSE,
              help = "validate inputs, print planned steps, write nothing")
)
opt <- parse_args(OptionParser(option_list = opt_list,
  description = "pi-matched enrichment of introgression outliers in selective sweeps"))

# ---- logging
logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
warnf <- function(...) { msg <- sprintf(...); logf("WARNING: %s", msg); warning(msg, call. = FALSE) }
log_filter <- function(what, n_in, n_out) logf("  filter %-45s in %7d  kept %7d  dropped %7d",
                                               what, n_in, n_out, n_in - n_out)
na_report <- function(what, x) {                  # every join reports its NAs
  n_na <- sum(is.na(x)); logf("  NA check %-43s %7d / %d NA", what, n_na, length(x)); invisible(n_na)
}

need <- function(x, flag) if (is.null(x) || is.na(x) || x == "") stop("Missing required argument ", flag, call. = FALSE)
need(opt$lineage, "--lineage")
if (!opt$lineage %in% c("Occidental", "Oriental")) stop("--lineage must be Occidental or Oriental")
if (!opt$pi_mode %in% c("aggregate", "bed"))      stop("--pi_mode must be aggregate or bed")
need(opt$dinv, "--dinv")
if (!opt$export_windows_only) {
  for (f in c("sweep_bed", "mask_bed", "genes", "pixy_pi")) need(opt[[f]], paste0("--", f))
}
TAG     <- opt$tag %||% paste(opt$lineage, opt$pi_mode, sep = ".")
out_path <- function(suffix) file.path(opt$out_dir, paste0(TAG, ".", suffix))

expand_glob <- function(x, what) {
  pats  <- str_trim(str_split(x, ",")[[1]])
  files <- unique(unlist(lapply(pats, Sys.glob)))
  if (length(files) == 0) stop("No file matched ", what, ": ", x, call. = FALSE)
  files
}
check_file <- function(f, what) if (!file.exists(f)) stop(what, " not found: ", f, call. = FALSE)

logf("s08.x1 pi-matched enrichment | lineage %s | pi_mode %s | tag %s", opt$lineage, opt$pi_mode, TAG)
logf("Recombination rate is not available: NOT included as a covariate (see header).")

# ================================================================
# 0. validate inputs
# ================================================================
dinv_files <- expand_glob(opt$dinv, "--dinv")
dinv_trio  <- str_match(basename(dinv_files), "([^.]+)_localFstats[^/]*\\.txt$")[, 2]
if (anyNA(dinv_trio)) stop("Cannot parse trio from: ", paste(dinv_files[is.na(dinv_trio)], collapse = ", "))
trios <- sort(unique(dinv_trio))
if (!is.null(opt$trios)) {
  keep  <- str_trim(str_split(opt$trios, ",")[[1]])
  miss  <- setdiff(keep, trios)
  if (length(miss)) stop("--trios not found among --dinv files: ", paste(miss, collapse = ", "))
  trios <- keep
}
trio_pops <- set_names(map(trios, ~ str_extract_all(.x, pattern_dict)[[1]]), trios)
bad <- names(which(lengths(trio_pops) != 3))
if (length(bad)) stop("Trio name does not resolve to 3 known pops: ", paste(bad, collapse = ", "))

detected <- if (any(grepl(OCC_PATTERN, trios))) "Occidental" else "Oriental"
if (detected != opt$lineage)
  stop("--lineage ", opt$lineage, " but trios look ", detected, " (rule of s07.x1): ",
       paste(trios, collapse = ", "), call. = FALSE)
logf("Trios (%d): %s", length(trios), paste(trios, collapse = ", "))

pop_map <- if (is.null(opt$pop_map)) {
  logf("No --pop_map: pixy pop names taken from the trio name.")
  tibble(trio = trios, P1 = map_chr(trio_pops, 1), P2 = map_chr(trio_pops, 2), P3 = map_chr(trio_pops, 3))
} else {
  check_file(opt$pop_map, "--pop_map")
  m <- read_tsv(opt$pop_map, col_types = cols(.default = "c"), comment = "#")
  if (!all(c("trio", "P1", "P2", "P3") %in% names(m))) stop("--pop_map needs columns trio, P1, P2, P3")
  miss <- setdiff(trios, m$trio)
  if (length(miss)) stop("--pop_map lacks trios: ", paste(miss, collapse = ", "))
  filter(m, trio %in% trios) %>% select(trio, P1, P2, P3)
}

if (!opt$export_windows_only) {
  check_file(opt$mask_bed, "--mask_bed"); check_file(opt$genes, "--genes")
  sweep_files <- expand_glob(opt$sweep_bed, "--sweep_bed")
  if (length(sweep_files) > 1) {                 # per-pop BEDs: keep this run's pops only (as s07.x1)
    sweep_pop <- map_chr(basename(sweep_files), ~ str_split(.x, "\\.")[[1]][3])
    keep <- sweep_pop %in% unique(unlist(trio_pops))
    logf("Sweep BEDs: %d of %d match the run's pops", sum(keep), length(keep))
    sweep_files <- sweep_files[keep]; sweep_pop <- sweep_pop[keep]
    if (!length(sweep_files)) stop("No sweep BED matched the run's pops.")
  }
  pixy_files  <- expand_glob(opt$pixy_pi, "--pixy_pi")
  if (!is.null(opt$chrom_sizes)) check_file(opt$chrom_sizes, "--chrom_sizes")
  if (!is.null(opt$s06_summary_dir) && !dir.exists(opt$s06_summary_dir))
    stop("--s06_summary_dir not found: ", opt$s06_summary_dir)
  if (!requireNamespace("future.apply", quietly = TRUE))
    stop("Package future.apply is required (install.packages('future.apply')).", call. = FALSE)
}

# ---- readers
read_bed3 <- function(f) {      # BED 0-based half-open -> 1-based GRanges
  d <- read_tsv(f, col_names = FALSE, col_types = cols(.default = "c"), comment = "#", progress = FALSE)
  if (ncol(d) < 3) stop("BED with < 3 columns: ", f)
  GRanges(d[[1]], IRanges(as.integer(d[[2]]) + 1L, as.integer(d[[3]])))
}
read_genes <- function(f) {
  if (grepl("\\.gff3?(\\.gz)?$", f, ignore.case = TRUE)) {
    g <- read_tsv(f, col_names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
                  col_types = "cccii----", comment = "#", progress = FALSE) %>% filter(type == "gene")
    logf("  genes: %d GFF 'gene' features from %s", nrow(g), basename(f))
    GRanges(g$chr, IRanges(g$start, g$end))
  } else {
    gr <- read_bed3(f); logf("  genes: %d BED intervals from %s", length(gr), basename(f)); gr
  }
}
read_dinv <- function(files) {
  map_dfr(files, function(f) {
    d <- read_tsv(f, col_types = cols(chr = "c", windowStart = "i", windowEnd = "i", .default = "d"),
                  progress = FALSE)
    miss <- setdiff(c("chr", "windowStart", "windowEnd", "D", "f_d", "f_dM", "d_f"), names(d))
    if (length(miss)) stop(basename(f), " lacks columns: ", paste(miss, collapse = ", "))
    d
  })
}
read_pixy <- function(files) {
  p <- map_dfr(files, ~ read_tsv(.x, col_types = cols(pop = "c", chromosome = "c", .default = "d"),
                                 na = c("NA", "nan", ""), progress = FALSE))
  miss <- setdiff(c("pop", "chromosome", "window_pos_1", "window_pos_2", "avg_pi", "no_sites",
                    "count_diffs", "count_comparisons", "count_missing"), names(p))
  if (length(miss)) stop("pixy files lack columns: ", paste(miss, collapse = ", "))
  p
}

# chromosome-name concordance: fail loudly
check_chr <- function(what, input_chr, win_chr, require_all_windows = FALSE) {
  input_chr <- unique(input_chr); win_chr <- unique(win_chr)
  shared <- intersect(input_chr, win_chr)
  logf("  chr check %-20s %3d chr, %3d shared with windows (%d window chr)",
       what, length(input_chr), length(shared), length(win_chr))
  if (length(shared) == 0 || length(shared) / length(input_chr) < 0.8)
    stop(what, " chromosome names do not match the Dinvestigate windows.\n  ", what, ": ",
         paste(head(input_chr, 5), collapse = ", "), "\n  windows: ",
         paste(head(win_chr, 5), collapse = ", "), call. = FALSE)
  if (require_all_windows && length(setdiff(win_chr, input_chr)))
    stop(what, " has no data for window chromosomes: ",
         paste(setdiff(win_chr, input_chr), collapse = ", "), call. = FALSE)
}

# ================================================================
# read windows (all trios); window_id is stable across runs
# ================================================================
logf("Reading Dinvestigate windows ...")
win <- map_dfr(trios, function(tr) {
  read_dinv(dinv_files[dinv_trio == tr]) %>% mutate(trio = tr, .before = 1)
}) %>%
  arrange(trio, chr, windowStart, windowEnd) %>%
  mutate(window_id = sprintf("%s|%s:%d-%d", trio, chr, windowStart, windowEnd), .before = 1)
if (anyDuplicated(win$window_id)) stop("Duplicated windows (same trio/chr/start/end) in --dinv input.")
walk(c("D", "f_d", "f_dM"), ~ na_report(paste("Dinvestigate", .x), win[[.x]]))
logf("  %d windows, %d chr", nrow(win), n_distinct(win$chr))
win_chr <- unique(win$chr)

# ---- dry run: validate the rest, print the plan
if (opt$dry_run) {
  if (!opt$export_windows_only) {
    pass_gr <- read_bed3(opt$mask_bed); check_chr("mask_bed", seqlevels(pass_gr), win_chr)
    for (f in sweep_files) check_chr(paste0("sweep:", str_split(basename(f), "\\.")[[1]][3]),
                                     seqlevels(read_bed3(f)), win_chr)
    check_chr("genes", seqlevels(read_genes(opt$genes)), win_chr)
    px <- read_pixy(pixy_files); check_chr("pixy", px$chromosome, win_chr, require_all_windows = TRUE)
    miss <- setdiff(unlist(pop_map[, c("P1", "P2", "P3")]), px$pop)
    if (length(miss)) stop("pixy pops missing: ", paste(miss, collapse = ", "))
    logf("  pixy: %d rows, pops %s, median window %s bp", nrow(px), paste(sort(unique(px$pop)), collapse = ","),
         format(median(px$window_pos_2 - px$window_pos_1 + 1)))
  }
  cat("\nPlanned steps (dry run, nothing written):\n",
      sprintf("  1 write %s and %s\n", out_path("dinv_windows.annot.bed"), out_path("dinv_windows.pixy_regions.1based.bed")),
      sprintf("  2 pi per window (P1/P2/P3), mode '%s'%s\n", opt$pi_mode,
              if (opt$pi_mode == "aggregate") sprintf(", drop windows with pixy coverage < %.2f", opt$min_pixy_cov) else ""),
      sprintf("  3 window table (PASS >= %.2f, in_sweep >= %d bp) -> %s\n", opt$pass_min_frac, opt$min_sweep_overlap, out_path("window_table.tsv.gz")),
      sprintf("  4 %d pi/length-matched controls per sweep window, caliper %.2f (log), %s -> %s\n",
              opt$n_match, opt$pi_tol, if (opt$no_thin) "no thinning" else "LD-thinned", out_path("matched_enrichment.tsv")),
      sprintf("  5 logistic regression per trio x test, %d circular sweep rotations, %d threads -> %s\n",
              opt$n_perm, opt$threads, out_path("regression.tsv")),
      sprintf("  6 window-length check -> %s\n", out_path("window_length_check.tsv")),
      sep = "")
  logf("Dry run OK.")
  quit(save = "no", status = 0)
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 1 export windows
# ================================================================
logf("STEP 1: export windows")
win %>%
  transmute(`#chr` = chr, start = windowStart - 1L, end = windowEnd, window_id, trio, D, f_d, f_dM) %>%
  write_tsv(out_path("dinv_windows.annot.bed"), na = "NA")
# pixy --bed_file: headerless, 3 columns, unique regions across trios. The pixy runs in
# run_Pixy_PiFstDxy* feed 1-based inclusive coordinates (*.pass.bed.1based.txt), so the
# pixy file keeps Dinvestigate's 1-based windowStart/windowEnd.
win %>% distinct(chr, windowStart, windowEnd) %>%
  write_tsv(out_path("dinv_windows.pixy_regions.1based.bed"), col_names = FALSE)
logf("  wrote %s (0-based BED, annotated) and %s (%d unique regions, 1-based, for pixy --bed_file)",
     basename(out_path("dinv_windows.annot.bed")), basename(out_path("dinv_windows.pixy_regions.1based.bed")),
     nrow(distinct(win, chr, windowStart, windowEnd)))
if (opt$export_windows_only) { logf("--export_windows_only: done."); quit(save = "no", status = 0) }

# ================================================================
# shared inputs: mask, sweeps, genes, chrom sizes
# ================================================================
logf("Reading mask / sweeps / genes ...")
pass_gr <- GenomicRanges::reduce(read_bed3(opt$mask_bed))
check_chr("mask_bed", seqlevels(pass_gr), win_chr)

gene_gr <- GenomicRanges::reduce(read_genes(opt$genes))
check_chr("genes", seqlevels(gene_gr), win_chr)

if (length(sweep_files) == 1) {
  logf("  one sweep BED: used as the target for every test (%s)", basename(sweep_files))
  single_sweep <- GenomicRanges::reduce(read_bed3(sweep_files))
  check_chr("sweep_bed", seqlevels(single_sweep), win_chr)
  single_sweep <- GenomicRanges::intersect(single_sweep, pass_gr)
  get_target <- function(pops) single_sweep
} else {
  sweep_list <- set_names(map(sweep_files, function(f) {
    gr <- GenomicRanges::reduce(read_bed3(f))
    check_chr(paste0("sweep:", str_split(basename(f), "\\.")[[1]][3]), seqlevels(gr), win_chr)
    GenomicRanges::intersect(gr, pass_gr)      # clip to PASS (idempotent if already clipped)
  }), sweep_pop)
  logf("  per-pop sweep BEDs: %s", paste(sprintf("%s(%d)", names(sweep_list), lengths(sweep_list)), collapse = " "))
  get_target <- function(pops) {               # shared sweep, as s07.x1
    pops <- intersect(pops, names(sweep_list))
    if (length(pops) < 2) return(GRanges())
    Reduce(GenomicRanges::intersect, sweep_list[pops])
  }
}

chr_len <- if (!is.null(opt$chrom_sizes)) {
  cs <- read_tsv(opt$chrom_sizes, col_names = FALSE, col_types = cols(.default = "c"), progress = FALSE)
  set_names(as.numeric(cs[[2]]), cs[[1]])
} else {
  logf("  no --chrom_sizes: chromosome length = max end over windows, mask, genes")
  bind_rows(tibble(chr = win$chr, end = win$windowEnd),
            tibble(chr = as.character(seqnames(pass_gr)), end = end(pass_gr)),
            tibble(chr = as.character(seqnames(gene_gr)), end = end(gene_gr))) %>%
    group_by(chr) %>% summarise(L = max(end), .groups = "drop") %>% deframe() %>% {set_names(as.numeric(.), names(.))}
}
if (length(setdiff(win_chr, names(chr_len)))) stop("--chrom_sizes lacks window chr: ",
                                                   paste(setdiff(win_chr, names(chr_len)), collapse = ", "))

# ---- GRanges helpers
cov_frac <- function(q_gr, s_gr) {            # fraction of each q range covered by (reduced) s_gr
  out <- numeric(length(q_gr))
  hits <- findOverlaps(q_gr, s_gr)
  if (length(hits)) {
    ov <- width(pintersect(q_gr[queryHits(hits)], s_gr[subjectHits(hits)]))
    s  <- rowsum(ov, queryHits(hits))
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out / width(q_gr)
}
thin_windows <- function(win_gr) {            # greedy non-overlapping set; same result as s07.x1,
  chr <- as.character(seqnames(win_gr))       # vectors pulled out of the GRanges for speed
  st  <- start(win_gr); en <- end(win_gr)
  ord <- order(chr, st)
  keep <- logical(length(win_gr)); last_end <- -Inf; last_chr <- ""
  for (i in ord) {
    if (chr[i] != last_chr || st[i] > last_end) {
      keep[i] <- TRUE; last_end <- en[i]; last_chr <- chr[i]
    }
  }
  keep
}
# rotate intervals along each chromosome by one offset per chr; lengths and spacing are
# kept, intervals crossing the chromosome end wrap to its start (split in two).
rotate_gr <- function(gr, offsets) {
  if (!length(gr)) return(gr)
  ch <- as.character(seqnames(gr))
  L  <- chr_len[ch]; k <- offsets[ch]
  s0 <- (start(gr) - 1 + k) %% L               # new 0-based start
  e0 <- s0 + width(gr)                         # new 0-based exclusive end
  wrap <- e0 > L
  c(GRanges(ch, IRanges(s0 + 1, pmin(e0, L))),
    GRanges(ch[wrap], IRanges(rep(1, sum(wrap)), e0[wrap] - L[wrap])))
}

# ================================================================
# 2  pi per window
# ================================================================
logf("STEP 2: pi per window (mode %s)", opt$pi_mode)
px <- read_pixy(pixy_files)
check_chr("pixy", px$chromosome, win_chr, require_all_windows = TRUE)
need_pops <- unique(unlist(pop_map[, c("P1", "P2", "P3")]))
miss <- setdiff(need_pops, px$pop)
if (length(miss)) stop("pixy output lacks pops: ", paste(miss, collapse = ", "), call. = FALSE)
px <- filter(px, pop %in% need_pops)
na_report("pixy count_comparisons", px$count_comparisons)

regions <- distinct(win, chr, windowStart, windowEnd)
reg_gr  <- GRanges(regions$chr, IRanges(regions$windowStart, regions$windowEnd))

pi_by_pop <- if (opt$pi_mode == "aggregate") {
  px_w  <- px$window_pos_2 - px$window_pos_1 + 1
  ratio <- median(px_w) / median(width(reg_gr))
  logf("  median pixy window %s bp, median Dinvestigate window %s bp (ratio %.1f)",
       format(median(px_w)), format(median(width(reg_gr))), ratio)
  if (ratio > opt$max_pixy_ratio) {
    msg <- sprintf(paste0("pixy windows are %.1fx longer than Dinvestigate windows: aggregated pi would be a ",
                          "regional, not a local, value. Re-run pixy on %s and use --pi_mode bed ",
                          "(or pixy --window_size <= %d)."),
                   ratio, basename(out_path("dinv_windows.pixy_regions.1based.bed")),
                   as.integer(opt$max_pixy_ratio * median(width(reg_gr))))
    if (!opt$allow_coarse_pi) stop(msg, " Override with --allow_coarse_pi.", call. = FALSE)
    warnf("%s", msg)
  }
  map_dfr(need_pops, function(p) {
    d  <- filter(px, pop == p, !is.na(count_comparisons), !is.na(count_diffs))
    logf("  pop %-12s pixy windows %d (dropped %d with NA counts)", p, nrow(d),
         sum(px$pop == p) - nrow(d))
    pgr <- GRanges(d$chromosome, IRanges(d$window_pos_1, d$window_pos_2))
    if (any(countOverlaps(pgr, pgr) > 1)) warnf("pixy windows overlap each other for pop %s", p)
    hits <- findOverlaps(reg_gr, pgr)
    qi <- queryHits(hits); si <- subjectHits(hits)
    ov <- width(pintersect(reg_gr[qi], pgr[si]))
    fr <- ov / width(pgr)[si]                  # fraction of the pixy window inside the overlap
    agg <- rowsum(cbind(d$count_diffs[si] * fr, d$count_comparisons[si] * fr, ov), qi)
    at  <- as.integer(rownames(agg))
    diffs <- comps <- cov_bp <- numeric(length(reg_gr))
    diffs[at] <- agg[, 1]; comps[at] <- agg[, 2]; cov_bp[at] <- agg[, 3]
    regions %>% mutate(pop = p,
                       pi = if_else(comps > 0, diffs / comps, NA_real_),   # never average avg_pi
                       pixy_cov = pmin(cov_bp / width(reg_gr), 1))
  })
} else {
  # pixy run on the exported regions: join on coordinates. pixy echoes the BED numbers,
  pxj <- px %>% transmute(pop, chr = chromosome, s = as.integer(window_pos_1), windowEnd = as.integer(window_pos_2),
                          pi = if_else(count_comparisons > 0, count_diffs / count_comparisons, NA_real_))
  n1 <- nrow(semi_join(pxj, regions, by = c("chr", "s" = "windowStart", "windowEnd")))
  n0 <- nrow(semi_join(mutate(pxj, s = s + 1L), regions, by = c("chr", "s" = "windowStart", "windowEnd")))
  shift <- if (n0 > n1) 1L else 0L
  logf("  pixy rows matching windows: %d as 1-based, %d as 0-based -> using %s", n1, n0,
       if (shift) "0-based" else "1-based")
  if (max(n0, n1) == 0) stop("No pixy row matches any window coordinate; was pixy run on the exported BED?")
  pxj <- pxj %>% mutate(windowStart = s + shift) %>% select(-s)
  dup <- pxj %>% count(pop, chr, windowStart, windowEnd) %>% filter(n > 1)
  if (nrow(dup)) stop(nrow(dup), " duplicated pixy rows per pop/window in bed mode.")
  map_dfr(need_pops, function(p)
    regions %>% left_join(filter(pxj, pop == p) %>% select(-pop), by = c("chr", "windowStart", "windowEnd")) %>%
      mutate(pop = p, pixy_cov = if_else(is.na(pi), 0, 1)))
}
pi_by_pop %>% group_by(pop) %>% group_walk(~ {
  logf("  pop %-12s windows %d | pi NA %d | pi == 0 %d | pixy_cov < %.2f: %d", .y$pop, nrow(.x),
       sum(is.na(.x$pi)), sum(.x$pi == 0, na.rm = TRUE), opt$min_pixy_cov, sum(.x$pixy_cov < opt$min_pixy_cov))
})

# ================================================================
# 3  analysis table
# ================================================================
logf("STEP 3: window table")
pi_pseudo <- if (is.na(opt$pi_pseudo)) min(pi_by_pop$pi[pi_by_pop$pi > 0], na.rm = TRUE) / 2 else opt$pi_pseudo
logf("  log(pi + %.3g) used for pi (pseudocount only matters for pi == 0 windows)", pi_pseudo)

TESTS <- tibble(test   = c("f_d_top", "f_dM_top", "f_dM_bottom"),
                flag   = c("out_fd_top", "out_fdm_top", "out_fdm_bot"),
                target = c("P2P3", "P2P3", "P1P3"),
                valid  = c("valid_fd", "valid_fdm", "valid_fdm"))

build_trio <- function(tr) {
  pm <- filter(pop_map, trio == tr)
  tp <- trio_pops[[tr]]
  d  <- filter(win, trio == tr)
  n0 <- nrow(d)

  # outlier flags: identical to s06 (quantiles over ALL windows of the trio), read by s07.x1
  valid_fd   <- !is.na(d$f_d) & d$f_d <= 1
  fd_clamped <- pmax(d$f_d, 0)
  q_fd <- if (any(valid_fd)) unname(quantile(fd_clamped[valid_fd], FD_TOP_Q)) else NA_real_
  valid_fdm <- !is.na(d$f_dM) & d$f_dM >= -1 & d$f_dM <= 1
  q_lo <- if (any(valid_fdm)) unname(quantile(d$f_dM[valid_fdm], 1 - FDM_TOP_Q)) else NA_real_
  q_hi <- if (any(valid_fdm)) unname(quantile(d$f_dM[valid_fdm], FDM_TOP_Q))     else NA_real_
  d <- d %>% mutate(valid_fd = valid_fd, valid_fdm = valid_fdm,
                    out_fd_top  = valid_fd  & !is.na(q_fd) & fd_clamped >= q_fd,
                    out_fdm_top = valid_fdm & !is.na(q_hi) & f_dM >= q_hi,
                    out_fdm_bot = valid_fdm & !is.na(q_lo) & f_dM <= q_lo)

  for (k in 1:3) {                             # pi of P1/P2/P3, NA-counted join
    pk <- pi_by_pop %>% filter(pop == pm[[paste0("P", k)]]) %>%
      select(chr, windowStart, windowEnd, pi, pixy_cov) %>%
      rename_with(~ paste0(.x, "_P", k), c(pi, pixy_cov))
    d <- left_join(d, pk, by = c("chr", "windowStart", "windowEnd"))
    if (nrow(d) != n0) stop("join on pi changed the row count for ", tr)
    na_report(sprintf("%s pi_P%d (%s)", tr, k, pm[[paste0("P", k)]]), d[[paste0("pi_P", k)]])
  }

  gr <- GRanges(d$chr, IRanges(d$windowStart, d$windowEnd))
  tgt <- list(P2P3 = get_target(c(tp[2], tp[3])), P1P3 = get_target(c(tp[1], tp[3])))
  d <- d %>% mutate(
    start = windowStart - 1L, end = windowEnd,              # BED coordinates
    window_length = windowEnd - windowStart + 1L,
    snp_density   = SNPS_PER_WINDOW / window_length,         # = 50/length: collinear, never both in a model
    gene_density  = cov_frac(gr, gene_gr),
    pass_fraction = cov_frac(gr, pass_gr),
    in_sweep_P2P3 = overlapsAny(gr, tgt$P2P3, minoverlap = opt$min_sweep_overlap),
    sweep_overlap_frac_P2P3 = cov_frac(gr, tgt$P2P3),
    in_sweep_P1P3 = overlapsAny(gr, tgt$P1P3, minoverlap = opt$min_sweep_overlap),
    sweep_overlap_frac_P1P3 = cov_frac(gr, tgt$P1P3),
    P1 = tp[1], P2 = tp[2], P3 = tp[3])

  # universe flags (no row is removed from the table; each filter is logged)
  ok_pass <- d$pass_fraction >= opt$pass_min_frac
  ok_cov  <- d$pixy_cov_P1 >= opt$min_pixy_cov & d$pixy_cov_P2 >= opt$min_pixy_cov & d$pixy_cov_P3 >= opt$min_pixy_cov
  ok_cov[is.na(ok_cov)] <- FALSE
  ok_pi   <- !is.na(d$pi_P1) & !is.na(d$pi_P2) & !is.na(d$pi_P3)
  logf("  trio %s  targets: P2P3 %d intervals (%.2f Mb), P1P3 %d intervals (%.2f Mb)", tr,
       length(tgt$P2P3), sum(width(tgt$P2P3)) / 1e6, length(tgt$P1P3), sum(width(tgt$P1P3)) / 1e6)
  log_filter(sprintf("%s PASS >= %.2f", tr, opt$pass_min_frac), n0, sum(ok_pass))
  log_filter(sprintf("%s + pixy coverage >= %.2f (P1,P2,P3)", tr, opt$min_pixy_cov), sum(ok_pass), sum(ok_pass & ok_cov))
  if (sum(ok_pass & !ok_cov)) warnf("%s: %d PASS windows dropped for pixy coverage < %.2f",
                                    tr, sum(ok_pass & !ok_cov), opt$min_pixy_cov)
  log_filter(sprintf("%s + pi not NA (P1,P2,P3)", tr), sum(ok_pass & ok_cov), sum(ok_pass & ok_cov & ok_pi))
  base <- ok_pass & ok_cov & ok_pi
  log_filter(sprintf("%s + f_d valid  -> universe f_d", tr), sum(base), sum(base & d$valid_fd))
  log_filter(sprintf("%s + f_dM valid -> universe f_dM", tr), sum(base), sum(base & d$valid_fdm))

  d %>% mutate(in_universe_fd = base & valid_fd, in_universe_fdm = base & valid_fdm) %>%
    select(window_id, trio, P1, P2, P3, chr, start, end, window_length, snp_density, gene_density,
           pass_fraction, pi_P1, pi_P2, pi_P3, pixy_cov_P1, pixy_cov_P2, pixy_cov_P3,
           in_sweep_P2P3, sweep_overlap_frac_P2P3, in_sweep_P1P3, sweep_overlap_frac_P1P3,
           D, f_d, f_dM, d_f, valid_fd, valid_fdm, out_fd_top, out_fdm_top, out_fdm_bot,
           in_universe_fd, in_universe_fdm)
}
tab <- map_dfr(trios, build_trio)
write_tsv(tab, out_path("window_table.tsv.gz"), na = "NA")
logf("  wrote %s (%d rows)", basename(out_path("window_table.tsv.gz")), nrow(tab))

# optional: recomputed flags vs s06 summary (s07.x1 reads the s06 flags)
if (!is.null(opt$s06_summary_dir)) {
  s06_cols <- c(out_fd_top  = paste0("If_top_",    100 - FD_TOP_Q  * 100, "_percent_f_d_outlier"),
                out_fdm_top = paste0("If_top_",    100 - FDM_TOP_Q * 100, "_percent_f_dM_outlier"),
                out_fdm_bot = paste0("If_bottom_", 100 - FDM_TOP_Q * 100, "_percent_f_dM_outlier"))
  walk(trios, function(tr) {
    f <- file.path(opt$s06_summary_dir, paste0(tr, "_LocalFstats_summary.txt"))
    if (!file.exists(f)) { logf("  s06 check: %s missing, skipped", basename(f)); return() }
    s <- read_tsv(f, show_col_types = FALSE, guess_max = Inf) %>%
      transmute(chr, start = windowStart - 1L, end = windowEnd, !!!syms(s06_cols))
    j <- inner_join(filter(tab, trio == tr), s, by = c("chr", "start", "end"), suffix = c("", ".s06"))
    mism <- map_int(names(s06_cols), ~ sum(j[[.x]] != j[[paste0(.x, ".s06")]], na.rm = TRUE))
    logf("  s06 check %s: %d/%d windows joined; flag mismatches fd_top %d, fdm_top %d, fdm_bot %d",
         tr, nrow(j), sum(tab$trio == tr), mism[1], mism[2], mism[3])
    if (any(mism > 0)) warnf("%s: recomputed outlier flags differ from s06", tr)
  })
}

# analysis frame for one trio x test (universe rows only)
test_frame <- function(tr, te) {
  u <- if (te$valid == "valid_fd") "in_universe_fd" else "in_universe_fdm"
  tab %>% filter(trio == tr, .data[[u]]) %>%
    transmute(window_id, chr, start, end,
              outlier  = .data[[te$flag]],
              in_sweep = .data[[paste0("in_sweep_", te$target)]],
              log_pi_P1 = log(pi_P1 + pi_pseudo), log_pi_P2 = log(pi_P2 + pi_pseudo),
              log_pi_P3 = log(pi_P3 + pi_pseudo), log_window_length = log(window_length),
              gene_density, pass_fraction)
}
COVS <- c("log_pi_P1", "log_pi_P2", "log_pi_P3", "log_window_length", "gene_density", "pass_fraction")
MATCH_VARS <- c("log_pi_P1", "log_pi_P2", "log_pi_P3", "log_window_length")

grid <- expand_grid(trio = trios, i_test = seq_len(nrow(TESTS)))

# ================================================================
# 4  pi-matched controls
# ================================================================
logf("STEP 4: pi-matched controls (n_match %d, caliper %.2f on log scale, %s)",
     opt$n_match, opt$pi_tol, if (opt$no_thin) "no thinning" else "LD-thinned windows")

match_one <- function(d, seed) {
  set.seed(seed)
  tr_i <- which(d$in_sweep); pool <- which(!d$in_sweep)
  X  <- as.matrix(d[, MATCH_VARS])
  Xp <- X[pool, , drop = FALSE]
  avail <- rep(TRUE, length(pool))
  n_got <- integer(length(tr_i)); ctrl_of <- vector("list", length(tr_i))
  for (j in sample(seq_along(tr_i))) {          # random order; nearest-neighbour, no replacement
    x   <- X[tr_i[j], ]
    dev <- abs(sweep(Xp, 2, x))
    ok  <- avail & rowSums(dev <= opt$pi_tol) == ncol(X)
    if (!any(ok)) next
    cand <- which(ok)
    dist <- sqrt(rowSums(dev[cand, , drop = FALSE]^2))
    pick <- cand[order(dist)][seq_len(min(opt$n_match, length(cand)))]
    avail[pick] <- FALSE
    n_got[j] <- length(pick); ctrl_of[[j]] <- pool[pick]
  }
  list(treated = tr_i, n_got = n_got, controls = ctrl_of)
}

balance_rows <- function(d, idx_s, idx_c, stage) {
  map_dfr(COVS, function(v) {
    xs <- d[[v]][idx_s]; xc <- d[[v]][idx_c]
    smd <- (mean(xs) - mean(xc)) / sqrt((var(xs) + var(xc)) / 2)
    bind_rows(tibble(group = "sweep", x = list(xs)), tibble(group = "control", x = list(xc))) %>%
      mutate(covariate = v, stage = stage,
             n = map_int(x, length), mean = map_dbl(x, mean), sd = map_dbl(x, sd),
             q25 = map_dbl(x, ~ unname(quantile(.x, 0.25))), median = map_dbl(x, median),
             q75 = map_dbl(x, ~ unname(quantile(.x, 0.75))), smd = smd) %>% select(-x)
  })
}

match_res <- pmap(grid, function(trio, i_test) {
  te <- TESTS[i_test, ]
  d  <- test_frame(trio, te)
  n_all <- nrow(d)
  if (!opt$no_thin) {
    d <- d[thin_windows(GRanges(d$chr, IRanges(d$start + 1L, d$end))), ]
    log_filter(sprintf("%s %s LD-thinning", trio, te$test), n_all, nrow(d))
  }
  id <- tibble(trio = trio, test = te$test, target = te$target)
  if (sum(d$in_sweep) == 0 || sum(!d$in_sweep) == 0) {
    warnf("%s %s: no windows inside (or outside) the %s target; matching skipped", trio, te$test, te$target)
    return(list(res = id, bal = NULL))
  }
  m <- match_one(d, opt$seed + i_test)
  matched <- m$n_got > 0
  s_idx <- m$treated[matched]; c_idx <- unlist(m$controls[matched])
  w_c   <- rep(1 / m$n_got[matched], m$n_got[matched])   # each sweep window weighs 1
  logf("  %s %-11s sweep windows %d: full %d, partial %d, unmatched %d; controls %d",
       trio, te$test, length(m$treated), sum(m$n_got == opt$n_match),
       sum(m$n_got > 0 & m$n_got < opt$n_match), sum(!matched), length(c_idx))
  if (sum(!matched)) warnf("%s %s: %d sweep windows had no control within the caliper",
                           trio, te$test, sum(!matched))

  a <- sum(d$outlier[s_idx]); b <- length(s_idx) - a
  cc <- sum(d$outlier[c_idx]); e <- length(c_idx) - cc
  mat <- matrix(c(a, b, cc, e), nrow = 2, byrow = TRUE)   # rows sweep/control, cols outlier yes/no
  fg <- tryCatch(fisher.test(mat, alternative = "greater"), error = function(err) NULL)
  f2 <- tryCatch(fisher.test(mat), error = function(err) NULL)

  res <- bind_cols(id, tibble(
    n_windows = nrow(d), n_sweep = length(m$treated),
    n_matched_full = sum(m$n_got == opt$n_match), n_matched_partial = sum(m$n_got > 0 & m$n_got < opt$n_match),
    n_unmatched = sum(!matched), n_controls = length(c_idx),
    sweep_outlier = a, sweep_nonoutlier = b, control_outlier = cc, control_nonoutlier = e,
    rate_sweep = a / (a + b), rate_control = cc / (cc + e),
    rate_control_weighted = sum(w_c * d$outlier[c_idx]) / sum(w_c),
    rate_sweep_unmatched_all = mean(d$outlier[m$treated]), rate_nonsweep_all = mean(d$outlier[!d$in_sweep]),
    odds_ratio = if (is.null(f2)) NA_real_ else unname(f2$estimate),
    or_ci_low  = if (is.null(f2)) NA_real_ else f2$conf.int[1],
    or_ci_high = if (is.null(f2)) NA_real_ else f2$conf.int[2],
    p_fisher_greater   = if (is.null(fg)) NA_real_ else fg$p.value,
    p_fisher_two_sided = if (is.null(f2)) NA_real_ else f2$p.value))
  bal <- bind_rows(balance_rows(d, which(d$in_sweep), which(!d$in_sweep), "before"),
                   balance_rows(d, s_idx, c_idx, "after")) %>%
    mutate(trio = trio, test = te$test, .before = 1)
  list(res = res, bal = bal)
})
matched_tbl <- map_dfr(match_res, "res") %>%
  mutate(p_fisher_greater_BH = p.adjust(p_fisher_greater, method = "BH"), .after = p_fisher_greater)
write_tsv(matched_tbl, out_path("matched_enrichment.tsv"), na = "NA")
write_tsv(map_dfr(match_res, "bal"), out_path("covariate_balance.tsv"), na = "NA")
logf("  wrote %s, %s", basename(out_path("matched_enrichment.tsv")), basename(out_path("covariate_balance.tsv")))

# ================================================================
# 5  logistic regression + circular permutation of sweeps
# ================================================================
# full model:  outlier ~ in_sweep + log(pi_P1) + log(pi_P2) + log(pi_P3) + log(window_length) + gene_density
# snp_density (= 50/window_length) is left out because it is collinear with log(window_length).

logf("STEP 5: logistic regression; in_sweep tested by U = O - E, %d circular permutations, %d threads",
     opt$n_perm, opt$threads)
logf("  minimum attainable permutation p = 1/(n_perm+1) = %.2g", 1 / (opt$n_perm + 1))
if (opt$threads > 1 && future::supportsMulticore()) {   # plan() substitutes its argument:
  future::plan("multicore", workers = opt$threads)       # keep the strategy a literal
} else if (opt$threads > 1) {
  future::plan("multisession", workers = opt$threads)
} else future::plan("sequential")

REG_VARS <- c("log_pi_P1", "log_pi_P2", "log_pi_P3", "log_window_length", "gene_density")
null_store <- list()

reg_tbl <- pmap_dfr(grid, function(trio, i_test) {
  te <- TESTS[i_test, ]
  tp <- trio_pops[[trio]]
  d  <- test_frame(trio, te)
  id <- tibble(trio = trio, test = te$test, target = te$target, n = nrow(d),
               n_outlier = sum(d$outlier), n_in_sweep = sum(d$in_sweep))
  if (sum(d$in_sweep) == 0 || sum(d$outlier) == 0) {
    warnf("%s %s: no in-sweep or no outlier windows; regression skipped", trio, te$test)
    return(id)
  }
  # null model -> expected outlier probability per window
  fit0 <- suppressWarnings(glm(reformulate(REG_VARS, "outlier"), data = d, family = binomial()))
  if (!fit0$converged) warnf("%s %s: covariate-only glm did not converge", trio, te$test)
  y  <- as.numeric(d$outlier)
  r0 <- y - fitted(fit0)                         # residual: observed - expected, per window
  O  <- sum(y[d$in_sweep]); E <- sum(fitted(fit0)[d$in_sweep]); U <- O - E
  separation <- O == 0 || O == sum(d$in_sweep)

  # full model: descriptive beta / Wald p (beta diverges under separation)
  fit <- suppressWarnings(glm(reformulate(c("in_sweep", REG_VARS), "outlier"), data = d, family = binomial()))
  cf  <- summary(fit)$coefficients
  b_obs <- cf["in_sweepTRUE", "Estimate"]
  if (separation) logf("  %s %s: %d/%d outliers in sweep -> complete separation, beta/Wald p not interpretable",
                       trio, te$test, O, sum(d$in_sweep))

  win_gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end))
  tgt <- get_target(if (te$target == "P2P3") c(tp[2], tp[3]) else c(tp[1], tp[3]))
  tgt_chr <- unique(as.character(seqnames(tgt)))
  chunks <- split(seq_len(opt$n_perm), cut(seq_len(opt$n_perm), max(1, opt$threads * 4), labels = FALSE))
  res <- future.apply::future_lapply(chunks, function(idx) {
    t(vapply(idx, function(b) {
      off <- set_names(floor(runif(length(tgt_chr)) * chr_len[tgt_chr]), tgt_chr)
      s   <- overlapsAny(win_gr, rotate_gr(tgt, off), minoverlap = opt$min_sweep_overlap)
      c(sum(r0[s]), sum(s))
    }, numeric(2)))
  }, future.seed = opt$seed + i_test)
  res  <- do.call(rbind, res)
  nv   <- res[, 1]
  n_empty <- sum(res[, 2] == 0)
  if (n_empty) logf("  %s %s: %d/%d rotations hit no window (U = 0, kept)", trio, te$test, n_empty, opt$n_perm)
  null_store[[paste(trio, te$test)]] <<- tibble(trio = trio, test = te$test, U = nv, observed = U)
  n_ge <- sum(nv >= U - 1e-9); n_le <- sum(nv <= U + 1e-9)   # tolerance: U is a sum of doubles
  logf("  %s %-11s O %d  E %.2f  O/E %.2f  U %.2f | beta %.3f (Wald p %.2g) | p_perm %.3g",
       trio, te$test, O, E, O / E, U, b_obs, cf["in_sweepTRUE", "Pr(>|z|)"], (n_ge + 1) / (opt$n_perm + 1))

  other <- as_tibble(as.list(set_names(cf[REG_VARS, "Estimate"], paste0("coef_", REG_VARS))))
  bind_cols(id, tibble(
    observed_in_sweep = O, expected_in_sweep = E, obs_over_exp = O / E, U_obs_minus_exp = U,
    p_perm_greater   = (n_ge + 1) / (opt$n_perm + 1),
    # doubled smaller tail (null of U is discrete and skewed for small sweeps)
    p_perm_two_sided = min(1, 2 * (min(n_ge, n_le) + 1) / (opt$n_perm + 1)),
    n_perm = opt$n_perm, n_perm_no_window = n_empty, p_perm_min = 1 / (opt$n_perm + 1),
    null_U_median = median(nv), null_U_mean = mean(nv), null_U_sd = sd(nv),
    separation = separation,
    coef_in_sweep = b_obs, se_in_sweep = cf["in_sweepTRUE", "Std. Error"],
    odds_ratio_in_sweep = exp(b_obs), p_wald = cf["in_sweepTRUE", "Pr(>|z|)"],
    converged = fit$converged), other)
})
reg_tbl <- reg_tbl %>%
  mutate(p_perm_greater_BH = p.adjust(p_perm_greater, method = "BH"), .after = p_perm_greater)
write_tsv(reg_tbl, out_path("regression.tsv"), na = "NA")
null_tbl <- bind_rows(null_store)
write_tsv(null_tbl, out_path("perm_null.tsv.gz"))
logf("  wrote %s, %s", basename(out_path("regression.tsv")), basename(out_path("perm_null.tsv.gz")))
future::plan("sequential")

# ================================================================
# 6  window length check
# ================================================================
logf("STEP 6: window length check")
logf("  NOTE: snp_density = %d / window_length, i.e. collinear with window_length; only log(window_length) is modelled.",
     SNPS_PER_WINDOW)
len_tbl <- expand_grid(trio = trios, target = c("P2P3", "P1P3")) %>% pmap_dfr(function(trio, target) {
  d <- tab %>% filter(.data$trio == !!trio, in_universe_fdm | in_universe_fd)
  s <- d[[paste0("in_sweep_", target)]]
  w <- if (any(s) && any(!s)) wilcox.test(d$window_length[s], d$window_length[!s]) else NULL
  cors <- map(1:3, ~ suppressWarnings(cor.test(d$window_length, d[[paste0("pi_P", .x)]], method = "spearman", exact = FALSE)))
  r <- tibble(trio = trio, target = target, n = nrow(d), n_in_sweep = sum(s),
              median_len_in = if (any(s)) median(d$window_length[s]) else NA_real_,
              median_len_out = median(d$window_length[!s]),
              p_wilcox = if (is.null(w)) NA_real_ else w$p.value,
              rho_len_pi_P1 = unname(cors[[1]]$estimate), rho_len_pi_P2 = unname(cors[[2]]$estimate),
              rho_len_pi_P3 = unname(cors[[3]]$estimate),
              p_rho_P1 = cors[[1]]$p.value, p_rho_P2 = cors[[2]]$p.value, p_rho_P3 = cors[[3]]$p.value)
  logf("  %s %s window_length median in %s / out %s bp (Wilcoxon p %.2g); Spearman rho(len, pi) P1 %.2f P2 %.2f P3 %.2f",
       trio, target, format(r$median_len_in), format(r$median_len_out), r$p_wilcox,
       r$rho_len_pi_P1, r$rho_len_pi_P2, r$rho_len_pi_P3)
  r
})
write_tsv(len_tbl, out_path("window_length_check.tsv"), na = "NA")

# =============
# plots
# ==============
logf("Plots")
theme_s07 <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "top")
plot_long <- tab %>% filter(in_universe_fdm | in_universe_fd) %>%
  select(trio, P1, P2, P3, window_length, in_sweep_P2P3, pi_P1, pi_P2, pi_P3) %>%
  pivot_longer(starts_with("pi_P"), names_to = "pop", values_to = "pi") %>%
  mutate(pop = sub("pi_", "", pop),
         pop = sprintf("%s (%s)", pop, case_when(pop == "P1" ~ P1, pop == "P2" ~ P2, TRUE ~ P3)),
         in_sweep = if_else(in_sweep_P2P3, "in P2&P3 sweep", "outside"),
         log10_pi = log10(pi + pi_pseudo))
sweep_cols <- c("in P2&P3 sweep" = COL_IN, "outside" = COL_OUT)

pdf(out_path("pi_in_vs_out_sweep.pdf"), width = 7, height = 3.2)
for (tr in trios) print(
  ggplot(filter(plot_long, trio == tr), aes(log10_pi, fill = in_sweep, colour = in_sweep)) +
    geom_density(alpha = 0.25, linewidth = 0.5) + facet_wrap(~ pop, nrow = 1) +
    scale_fill_manual(values = sweep_cols, name = NULL) + scale_colour_manual(values = sweep_cols, name = NULL) +
    labs(x = expression(log[10]~pi), y = "density", title = tr) + theme_s07)
invisible(dev.off())

pdf(out_path("window_length_vs_pi.pdf"), width = 7, height = 3.2)
for (tr in trios) print(
  ggplot(filter(plot_long, trio == tr), aes(log10(window_length), log10_pi, colour = in_sweep)) +
    geom_point(data = ~ filter(.x, in_sweep == "outside"), size = 0.4, alpha = 0.3, stroke = 0) +
    geom_point(data = ~ filter(.x, in_sweep != "outside"), size = 1.2, alpha = 0.9, stroke = 0) +
    facet_wrap(~ pop, nrow = 1) +
    scale_colour_manual(values = sweep_cols, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = expression(log[10]~"window length (bp)"), y = expression(log[10]~pi), title = tr) + theme_s07)
invisible(dev.off())

pdf(out_path("in_sweep_U_null.pdf"), width = 7, height = 2.6 * max(1, length(trios)))
if (nrow(null_tbl)) print(
  ggplot(null_tbl, aes(U)) +
    geom_histogram(bins = 60, fill = COL_OUT, colour = NA) +
    geom_vline(data = distinct(null_tbl, trio, test, observed), aes(xintercept = observed),
               colour = COL_IN, linewidth = 0.8) +
    facet_grid(trio ~ test, scales = "free") +
    labs(x = "U = observed - expected outliers in sweeps, under circular rotation (line: observed)", y = "permutations") +
    theme_s07)
invisible(dev.off())

logf("Done. Outputs in %s with tag %s", opt$out_dir, TAG)
