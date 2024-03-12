rm(list = ls())
library(ggplot2)
library(cowplot)

data_all <- read.table("fastst_cv.data", header = T, sep = "\t")

data_summary <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE)
    )
  }
  data_sum <- ddply(data, groupnames,
    .fun = summary_func,
    varname
  )
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(data_all,
  varname = "CV_error",
  groupnames = c("K")
)



p <- ggplot(df2, aes(y = CV_error, x = K)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = CV_error - sd, ymax = CV_error + sd),
    width = .2,
    position = position_dodge(0.05)
  ) +
  scale_x_continuous(breaks = c(1:20)) +
  xlab("K ancestry") + ylab("Cross balidation error") +
  theme_half_open()

pdf(file = "fastst_cv.pdf")

dev.off()
