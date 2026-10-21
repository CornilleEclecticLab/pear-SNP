# Downsampling robustness check

# Run with: Rscript s10.x3.compare_arms.R

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

if (interactive()) {
    script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
    setwd(script_path)
    rm(list=ls())
}

out_dir <- "../output/s12.downsample_AD"
YVAR    <- "Deleterious_per_synonymous_SNP"

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

arm_files <- list.files(out_dir, pattern = "\\.synonymous\\.pop_mutation_burden\\.Rdata\\.txt$",
                        full.names = TRUE)
names(arm_files) <- sub("\\.synonymous.*", "", basename(arm_files))

df <- bind_rows(lapply(names(arm_files), function(a)
        read.table(arm_files[[a]], header = TRUE) %>% mutate(Arm = a))) %>%
  mutate(pop  = factor(Population, levels = desired_levels),
         side = if_else(pop %in% west_pops, "West", "East"),
         crop_type = case_when(pop %in% wild_pops ~ "Wild",
                               pop %in% crop_pops ~ "Cultivar",
                               TRUE ~ "other")) %>%
  filter(crop_type %in% c("Wild", "Cultivar")) %>%
  mutate(crop_type = relevel(factor(crop_type), ref = "Wild"))

test_one <- function(dd) {
  m0 <- lmer(reformulate(c("crop_type", "(1 | pop)"), YVAR), data = dd)
  cf <- coef(summary(m0))["crop_typeCultivar", ]
  data.frame(
    n_Cultivar      = sum(dd$crop_type == "Cultivar"),
    n_Wild          = sum(dd$crop_type == "Wild"),
    median_Cultivar = median(dd[[YVAR]][dd$crop_type == "Cultivar"]),
    median_Wild     = median(dd[[YVAR]][dd$crop_type == "Wild"]),
    wilcox_p        = wilcox.test(reformulate("crop_type", YVAR), data = dd)$p.value,
    lmm_estimate    = cf[["Estimate"]],
    lmm_p           = cf[["Pr(>|t|)"]]
  )
}

res <- df %>%
  group_by(Arm, side, Allelic_state) %>%
  group_modify(~ test_one(.x)) %>%
  ungroup() %>%
  mutate(direction = if_else(median_Cultivar < median_Wild, "Cultivar<Wild", "Cultivar>Wild"))

# Change of the cultivar effect relative to recall_full (same caller, full depth)
res <- res %>%
  group_by(side, Allelic_state) %>%
  mutate(lmm_change_vs_recall_full_pct =
           100 * (lmm_estimate - lmm_estimate[Arm == "recall_full"]) /
           abs(lmm_estimate[Arm == "recall_full"])) %>%
  ungroup() %>%
  arrange(side, Allelic_state, Arm)

write.table(res, file.path(out_dir, "s10.compare_arms.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
