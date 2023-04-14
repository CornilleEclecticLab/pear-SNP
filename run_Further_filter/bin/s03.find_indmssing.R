rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)

sumfile <- "../output/whole_pear/whole_pear.combine.sum_number.txt"
ind_miss_file <- "../output/whole_pear/whole_pear.combine.ind_miss.txt"

nu <- read.table(sumfile)

sum <- nu$V1

miss_df <- read.table(ind_miss_file)
names(miss_df) <- c("individual", "missing_number")

miss_df$missing_rate <- (miss_df$missing_number / sum) * 100


p1 <- ggplot(miss_df, aes(x = missing_rate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "yellow") +
  xlab("missiong rate (%)") + 
  theme_cowplot(12)

p3 <- ggplot(miss_df, aes(x = reorder(individual, missing_rate), y = missing_rate)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_line(yintercept = 30, linetype = "dashed", color = "yellow") +
  xlab("individual") +
  ylab("missiong rate (%)") +
  coord_flip() +
  theme_cowplot(12)

ggsave("missing_Hist.png", plot = p1, dpi = 300)
ggsave("missing_bar.png", plot = p3, dpi = 300, 
       width = 8,
       height = 45,
       units = "in")

p2 <- ggplot(miss_df, aes(x = missing_rate, y = "SNP")) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "yellow") +
  theme_cowplot(12)

remove_df <- filter(miss_df, miss_df$missing_rate > 30)

fileConn <- file("remove.list")
writeLines(as.character(remove_df$individual), fileConn)
close(fileConn)
