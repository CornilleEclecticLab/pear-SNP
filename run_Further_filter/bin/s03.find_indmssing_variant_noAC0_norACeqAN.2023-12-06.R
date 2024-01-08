rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)


setwd("/shared/home/ynie/work/pear/run_Further_filter/bin/")

sumfile <- "/shared/ifbstor1/projects/pear_snp3/pear/run_Further_filter/output/s04./pear_Oct2023.Combine_chr.vcf.gz.SNPs_number_txt"
ind_miss_file <- "/shared/ifbstor1/projects/pear_snp3/pear/run_Further_filter/output/s04./pear_Oct2023.Combine_chr.vcf.gz.ind_miss_txt"

nu <- read.table(sumfile)

sum <- nu$V1

miss_df <- read.table(ind_miss_file)
names(miss_df) <- c("individual", "missing_number","species","depth", "dataset")

miss_df$missing_rate <- (miss_df$missing_number / sum) * 100


# Filter species with fewer than 10 individuals
filtered_miss_df <- miss_df %>%
  group_by(species) %>%
  filter(n() >= 12)


# Filter species with fewer than 10 individuals
filtered_missing_rate_miss_df <- filtered_miss_df %>%
  filter(missing_rate <= 40)


p1 <- ggplot(miss_df, aes(x = missing_rate)) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "blue") +
  # geom_vline(xintercept = 30, linetype = "dashed", color = "green") +
  xlab("missing rate per sample (%)") + 
  theme_bw(base_size = 12)
  # theme_cowplot(12)
p1

# p3 <- ggplot(miss_df, aes(x = reorder(individual, missing_rate), y = missing_rate)) +
  p3 <- ggplot(miss_df, aes(x = individual, y = missing_rate, fill = species)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "green") +
  xlab("individual") +
  ylab("missing rate (%)") +
  coord_flip() +
  theme_bw(base_size = 12)
  # theme_cowplot(12)
p3


p3 <- ggplot(filtered_miss_df, aes(x = individual, y = depth, fill = dataset, group = interaction(dataset, species))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 40, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "green") +
  geom_text(aes(label = sprintf("%.2f",missing_rate)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +  # Add this line for text labels
  xlab("individual") +
  ylab("Depth(% mising rate in digital)") +
  coord_flip() +
  theme_bw(base_size = 12)

print(p3)



ggsave("missing_bar_depth.pdf", plot = p3, dpi = 300, 
       width = 8,
       height = 88,
       units = "in",
       limitsize = FALSE)




p4 <- ggplot(filtered_missing_rate_miss_df, aes(x = individual, y = missing_rate, fill = dataset, group=interaction(dataset, species))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "green") +
  geom_hline(yintercept = 40, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.1f",depth)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") + 
  ylab("Missing Rate (%)") +
  coord_flip() +
  facet_wrap(~ species, scales = "free_y") +  # Facet by species
  theme_bw(base_size = 12)

print(p4)

ggsave("missing_hist.overlap.missingRateLower30.pdf", plot = p4, dpi = 300,
       width = 13,
       height = 30,
       units = "in",
       limitsize = FALSE)




# Create violin plot for the filtered data
p5 <- ggplot(filtered_missing_rate_miss_df, aes(x = species, y = missing_rate, fill = species)) +
  geom_violin() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "green") +
  xlab("Species") +
  ylab("Missing Rate (%)") +
  theme_bw(base_size = 12)

print(p5)

ggsave("missing_Violin.png", plot = p5, dpi = 300)



ggsave("missing_Hist.png", plot = p1, dpi = 300)

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

