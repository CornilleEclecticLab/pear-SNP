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



data2 <- read.table("../input_local/depth.kinship_cutoff0.125.txt",
                   header = TRUE)
head(data2)

p <- ggplot() +
  geom_histogram(data=data2, aes(x=depth), bins =40,  
                 fill="lightcoral", color="black",alpha=0.5)+
  ylab("Num of pear individuals from Zhang2021" )+
  xlab("Depth  (times based on 540M bp)")+
  geom_vline(aes(xintercept=10))+
  theme_classic()+
  scale_x_continuous(breaks=seq(5,16,1)) 
p

p+geom_density(data=data2, aes(x=depth)) 

p



data3 <- read.table("../input_local/depth.project.removeUnclearSpeciesPrunus.RemoveGroupLessThan8individuals.txt",
                    header = TRUE)
p <- ggplot()+
  geom_histogram(data=data3, aes(x=Depth, colors="black",fill=Dataset))+
  theme_classic()+
  geom_density(alpha=0.4)
  scale_x_continuous(breaks=seq(0,40,5)) 
p

p <- ggplot()+
  geom_histogram(data=data3, aes(x=Depth, color=Group))+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,40,5)) 
p
