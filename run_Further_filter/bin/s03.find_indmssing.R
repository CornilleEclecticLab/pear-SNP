rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)

sumfile <- "/shared/home/ynie/work/pear/run_Further_filter/output/pear_Jul2023.Fmissing02_m2M2_maff005/pear_Jul2023.Fmiss02_m2M2_maf005.sum_number.txt"
ind_miss_file <- "/shared/home/ynie/work/pear/run_Further_filter/output/pear_Jul2023.Fmissing02_m2M2_maff005/pear_Jul2023.Fmiss02_m2M2_maf005.ind_miss.txt"

nu <- read.table(sumfile)

sum <- nu$V1

miss_df <- read.table(ind_miss_file)
names(miss_df) <- c("individual", "missing_number")

miss_df$missing_rate <- (miss_df$missing_number / sum) * 100


p1 <- ggplot(miss_df, aes(x = missing_rate)) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "yellow") +
  xlab("missing rate (%)") + 
  theme_bw(base_size = 12)
  # theme_cowplot(12)


p3 <- ggplot(miss_df, aes(x = reorder(individual, missing_rate), y = missing_rate)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  geom_line(yintercept = 30, linetype = "dashed", color = "yellow") +
  xlab("individual") +
  ylab("missing rate (%)") +
  coord_flip() +
  theme_bw(base_size = 12)
  # theme_cowplot(12)

ggsave("missing_Hist.png", plot = p1, dpi = 300)
ggsave("missing_bar.png", plot = p3, dpi = 300, 
       width = 8,
       height = 66,
       units = "in",
       limitsize = FALSE)

p2 <- ggplot(miss_df, aes(x = missing_rate, y = "SNP")) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "yellow") +
  theme_bw(base_size = 12)
  #theme_cowplot(12)

remove_df <- filter(miss_df, miss_df$missing_rate > 30)

fileConn <- file("remove.list")
writeLines(as.character(remove_df$individual), fileConn)
close(fileConn)

