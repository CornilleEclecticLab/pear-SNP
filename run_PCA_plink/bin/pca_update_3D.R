rm(list = ls())

library("ggplot2")
library("scatterplot3d")

data <- read.table("./K6M_pca_shape.txt", header = T)

pdf("PCA_shape.pdf", width=7, height=6)
p3d <- scatterplot3d(data[, c(11, 12, 13)],
  pch = data$shape,
  color = data$color,
  angle = 50,
  xlab = "PC1 (58.4164% explained)",
  ylab = "PC2 (50.6314% explained)",
  zlab = "PC3 (10.83% explained)"
)

legend("right",
  legend = c(
    "DOM_CID", "DOM_DES",
    "SIE_KAZ", "SIE_CHN", "ORI_ARM",
    "SYL_DNK", "SYL_FRA_WEST", "SYL_FRA_EAST", "SYL_DEU", "SYL_AUT", "SYL_ROU", "SYL_ITA", "SYL_MKD",
    "BAC_CHN"
  ),
  col = c(
    "#ffff44", "#ff942d",
    "#ff2920", "#ff2920", "#ff2920",
    "#002cf9", "#002cf9", "#002cf9", "#002cf9", "#8e8e93", "#9b1e8f", "#8e8e93", "#9b1e8f",
    "#00fb3a"
  ),
  pch = c(
    0, 1,
    2, 3, 4, 
    5,     6, 7, 8, 9, 10, 11, 12,
    13
  )
)
dev.off()



# legend("right",
#   legend = c(
#     "DOM_CID", "DOM_DES",
#     "SIE_KAZ", "SIE_CHN", "ORI_ARM",
#     "SYL_DNK", "SYL_FRA_WEST", "SYL_FRA_EAST", "SYL_DEU", "SYL_AUT", "SYL_ROU", "SYL_ITA", "SYL_MKD",
#     "BAC_CHN",
#     "Admix_Others"
#   ),
#   col = c(
#     "#ffff44", "#ff942d",
#     "#ff2920", "#ff2920", "#ff2920",
#     "#002cf9", "#002cf9", "#002cf9", "#002cf9", "#8e8e93", "#9b1e8f", "#8e8e93", "#9b1e8f",
#     "#00fb3a",
#     "#000000"
#   ),
#   pch = c(
#     15, 15,
#     15, 16, 17,
#     15, 16, 17, 18, 5, 15, 6, 16,
#     11,
#     4
#   )
# )


data <- read.table("./K6M_pca_SYL_shape.txt", header = T)

pdf("PCA_SYL_shape.pdf", width=7, height=6)
p3d <- scatterplot3d(data[, c(11, 12, 13)],
  pch = data$shape,
  color = data$color,
  xlim = c(-0.04, 0),
  angle = 50,
  xlab = "PC1 (58.4164% explained)",
  ylab = "PC2 (50.6314% explained)",
  zlab = "PC3 (10.83% explained)"
)

legend("right",
  legend = c(
    "SYL_DNK", "SYL_FRA_WEST", "SYL_FRA_EAST", "SYL_DEU", "SYL_AUT", "SYL_ROU", "SYL_ITA", "SYL_MKD"
  ),
  col = c(
    "#002cf9", "#002cf9", "#002cf9", "#002cf9", "#8e8e93", "#9b1e8f", "#9b1e8f", "#8e8e93"
  ),
  pch = c(
    5, 6, 7, 8, 9, 10, 11, 12
  )
)
dev.off()

# legend("right",
#   legend = c(
#     "SYL_DNK", "SYL_FRA_WEST", "SYL_FRA_EAST", "SYL_DEU", "SYL_AUT", "SYL_ROU", "SYL_MKD", "SYL_ITA",
#     "Admix_Others"
#   ),
#   col = c(
#     "#002cf9", "#002cf9", "#002cf9", "#002cf9", "#8e8e93", "#9b1e8f", "#9b1e8f", "#8e8e93",
#     "#000000"
#   ),
#   pch = c(
#     15, 16, 17, 18, 5, 15, 16, 6,
#     4
#   )
# )
