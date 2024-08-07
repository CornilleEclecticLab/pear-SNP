library(ggplot2)
library(ggrepel)
# install.packages("ggforce")
library(ggforce)
library(adegenet)

setwd("/shared/ifbstor1/projects/pear_snp3/pear/run_PCA_plink/output/branch15.batch24.pear_Dec2023.noClone")

rm(list=ls())
# Read the data file
data <- read.table("input_pca.data",sep ='\t', header = TRUE)

# Convert Group and Color columns to factors
data$Group <- as.factor(data$Group)
data$Color <- as.factor(data$Color)

# Set an option for avoiding label overlaps
options(ggrepel.max.overlaps = 300)

# Create the ggplot object


p02 <- ggplot(data, aes(x = PC1, y = PC2,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1 (29.8%) explained") +
  ylab("PC2 (7.7%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()
p02



p03 <- ggplot(data, aes(x = PC3, y = PC4,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC3 (5.0%) explained") +
  ylab("PC4 (3.0%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()
p03


# Plot with lable and segments
# Function to calculate central points for each group
calculate_central_points <- function(group_data) {
  central_point <- colMeans(group_data[, c("PC1", "PC2")])
  return(data.frame(Group = unique(group_data$Group), Central_PC1 = central_point[1], Central_PC2 = central_point[2], Color_c = group_data$Color[1]))
}

# Calculate central points for each group
central_points <- do.call(rbind, by(data, data$Group, calculate_central_points))

# Merge central points back into the PCA data
data <- merge(data, central_points, by = "Group")


p04 <- ggplot(data, aes(x = PC1, y = PC2,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1 (29.8%) explained") +
  ylab("PC2 (7.7%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()+
  geom_segment(aes(x = Central_PC1, y = Central_PC2, xend = PC1, yend = PC2), alpha = 0.5) +
  geom_label(aes(x = Central_PC1, y = Central_PC2, label = Group), alpha = 0.2)
p04

p05 <- ggplot(data, aes(x = PC1, y = PC2,color=Color, shape = Group,label = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1 (29.8%) explained") +
  ylab("PC2 (7.7%) explained") +
  theme_classic() +
  scale_color_identity()+
  coord_fixed()+
  geom_segment(aes(x = Central_PC1, y = Central_PC2, xend = PC1, yend = PC2), alpha = 0.5)

p05 <- p05 + geom_label(data = central_points, aes(x = Central_PC1, y = Central_PC2, label = Group, color=Color_c), alpha = 1)

print(p05)




pca_result <- prcomp(data[c(4:15)])
summary(pca_result)$importance

scaled_data <- scale(data[c(4:15)])


# Following code is from mapmixture, the function of calculating centroid is the same with before

# Calculate centroid (mean average) positions for each group
centroid_df <- stats::aggregate(
  cbind(PC1, PC2) ~ Group,
  FUN = mean,
  data = data
)

# Edit column names of centroid data frame
colnames(centroid_df)[2:3] <- c("cen1", "cen2")

# Add centroid coordinates to main data frame
df <- dplyr::left_join(data, centroid_df, by = "Group")

