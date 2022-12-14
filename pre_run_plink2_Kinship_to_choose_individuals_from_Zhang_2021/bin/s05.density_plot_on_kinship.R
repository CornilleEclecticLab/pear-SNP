rm(list = ls())
library(ggplot2)


data <- read.table("../output/without_clone/Pear.cutoff0.354.kin0",
                   header = FALSE)

head(data)
kinship <- data.frame(data$V8)

p <- ggplot() +
  geom_density(data=data, aes(x=V8)) +
  scale_x_continuous(breaks=seq(-1,0.5,0.1)) +
  ylab("Density of pear from Zhang2021" )+
  xlab("Kinship")+
  geom_vline(aes(xintercept=0.354))
p

