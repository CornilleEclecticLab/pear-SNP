# Depth control for deleterious burden
#
# 
# Mixed models with sequencing depth as a fixed covariate and population as a
# random intercept:
#   m0   : burden ~ crop_type + (1|pop)                          (no depth, baseline)
#   m1   : burden ~ crop_type + depth + (1|pop)                  (depth covariate)
#   m2   : burden ~ crop_type + depth_within + (1|pop)           (within-population depth effect)
# P-values from lmerTest (Satterthwaite df). Run with R 4.3.1.

# Depth = mean mapped depth after MarkDuplicates

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}

out_dir <- "../output/s11.depth_control"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#############
# Populations
#############
desired_levels <- c(
  "pyra", "cauc",
  "comm_Dessert", "comm_Perry",
  "pash", "ussu", "betu",
  "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White"
)
west_pops <- desired_levels[1:4]
east_pops <- desired_levels[5:12]
wild_pops <- c("pyra", "cauc", "pash", "ussu")
crop_pops <- c("comm_Dessert", "comm_Perry",
               "Sand_CN", "Sand_CN-SW", "Sand_CN-SE", "pyri_JP", "White")

###########
# Read data
###########
df_burden <- read.table("../output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt",
                        header = TRUE)

depth_ppy <- read.delim("depth.txt", check.names = FALSE) %>%
  select(Sample_ID = ID_vcf, Depth_PPY = `MapppedDetpth(X)`)

depth_pcom <- read.table("../../run_GATK_variant_calling_refer_on_communis/output/s02.stat_format.pear.ref_comm.Jul2024.txt",
                         header = TRUE) %>%
  select(Sample_ID = Sample, Depth_PCOM = Markdu_avg_depth)

df <- df_burden %>%
  left_join(depth_ppy,  by = "Sample_ID") %>%
  left_join(depth_pcom, by = "Sample_ID") %>%
  mutate(pop  = factor(Population, levels = desired_levels),
         side = case_when(pop %in% west_pops ~ "West",
                          pop %in% east_pops ~ "East"),
         crop_type = case_when(pop %in% wild_pops ~ "Wild",
                               pop %in% crop_pops ~ "Cultivar",
                               TRUE ~ "other"),
         Depth = if_else(side == "West", Depth_PCOM, Depth_PPY)) %>%
  filter(crop_type %in% c("Wild", "Cultivar"))   # rootstock (betu) excluded
stopifnot(!any(is.na(df$Depth)))

##############################
# Depth summary per population
##############################
depth_by_pop <- df %>%
  distinct(Sample_ID, side, pop, crop_type, Depth) %>%
  group_by(side, crop_type, pop) %>%
  summarise(n = n(), median = median(Depth), mean = mean(Depth),
            min = min(Depth), max = max(Depth), .groups = "drop") %>%
  arrange(side, pop)
print(depth_by_pop, n = Inf)
write.table(depth_by_pop, file.path(out_dir, "s11.depth_by_population.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

depth_test <- df %>%
  distinct(Sample_ID, side, crop_type, Depth) %>%
  group_by(side) %>%
  summarise(median_Cultivar = median(Depth[crop_type == "Cultivar"]),
            median_Wild     = median(Depth[crop_type == "Wild"]),
            wilcox_p        = wilcox.test(Depth ~ crop_type)$p.value,
            .groups = "drop")
print(depth_test)
write.table(depth_test, file.path(out_dir, "s11.depth_cultivar_vs_wild.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

########
# Models
########
tidy_lmer <- function(m, model_name) {
  cf <- coef(summary(m))
  data.frame(model = model_name, term = rownames(cf),
             estimate = cf[, "Estimate"], se = cf[, "Std. Error"],
             df = cf[, "df"], p = cf[, "Pr(>|t|)"],
             singular = isSingular(m), row.names = NULL)
}

fit_models <- function(dd) {
  dd <- dd %>%
    mutate(crop_type = relevel(factor(crop_type), ref = "Wild"),
           pop       = droplevels(pop),
           depth     = Depth - mean(Depth)) %>%
    group_by(pop) %>%
    mutate(depth_within = Depth - mean(Depth)) %>%
    ungroup()

  y <- "Deleterious_per_synonymous_SNP"
  m0 <- lmer(reformulate(c("crop_type", "(1 | pop)"), y), data = dd)
  m1 <- lmer(reformulate(c("crop_type", "depth", "(1 | pop)"), y), data = dd)
  m2 <- lmer(reformulate(c("crop_type", "depth_within", "(1 | pop)"), y), data = dd)

  bind_rows(tidy_lmer(m0, "m0_no_depth"),
            tidy_lmer(m1, "m1_depth"),
            tidy_lmer(m2, "m2_depth_within"))
}

res <- df %>%
  group_by(side, Allelic_state, Region_type) %>%
  group_modify(~ fit_models(.x)) %>%
  ungroup() %>%
  mutate(p_fmt = format.pval(p, digits = 3))

write.table(res, file.path(out_dir, "s11.mixed_model_all_terms.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################
# Summary: cultivar effect with/without depth, and depth effects
################################################################
key <- res %>%
  filter(term %in% c("crop_typeCultivar", "depth", "depth_within")) %>%
  mutate(col = paste(model, term, sep = ":")) %>%
  select(side, Allelic_state, Region_type, col, estimate, p) %>%
  pivot_wider(names_from = col, values_from = c(estimate, p), names_glue = "{col}.{.value}") %>%
  mutate(cultivar_effect_change_pct =
           100 * (`m1_depth:crop_typeCultivar.estimate` - `m0_no_depth:crop_typeCultivar.estimate`) /
           abs(`m0_no_depth:crop_typeCultivar.estimate`)) %>%
  arrange(side, Region_type, Allelic_state)

write.table(key, file.path(out_dir, "s11.mixed_model_key_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(key %>%
        filter(Region_type == "all_masked") %>%
        select(side, Allelic_state,
               `m0_no_depth:crop_typeCultivar.estimate`, `m0_no_depth:crop_typeCultivar.p`,
               `m1_depth:crop_typeCultivar.estimate`, `m1_depth:crop_typeCultivar.p`,
               `m1_depth:depth.p`, `m2_depth_within:depth_within.p`,
               cultivar_effect_change_pct),
      width = Inf)
