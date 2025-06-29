library(ggplot2)
library(readr)
library(dplyr)


df <- read_tsv("fmiss.txt", col_names = c("CHROM", "POS", "F_MISSING"))


ggplot(df, aes(x = F_MISSING)) +
  #geom_histogram(binwidth = 0.01, fill = "#69b3a2", color = "white", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  labs(title = "F_MISSING distribution",
       x ="F_MISSING",
       y = "frequency") +
  theme_minimal()



ggplot(df, aes(x = F_MISSING)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001,
                 fill = "#69b3a2", color = "white", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  labs(title = "F_MISSING distribution",
       x ="F_MISSING",
       y = "frequency") +
  theme_minimal()
