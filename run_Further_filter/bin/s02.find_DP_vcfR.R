rm(list = ls())
library(vcfR)
library(reshape2)
library(ggplot2)

rand_vcf <- "../output/pear.test_rand.new_17scaffolds.vcf.gz"

vcf <- read.vcfR(rand_vcf) 

dp <- extract.gt(vcf, element="DP", as.numeric = T)

dpf <- melt(dp,varnames=c("Index","Sample"),value.name="Depth",na.rm=TRUE)
# dpf <- dpf[ dpf$Depth > 0,]


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

p <- ggplot(dpf, aes(x=Sample,y=Depth)) + 
  geom_violin(fill="#C0C0C0", adjust=1.0, scale="count",trim=TRUE) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust =1)) +
  guides(fill="none") +
  stat_summary(fun.data = quantiles_95, geom="boxplot", width=0.5, alpha=0.6) +
  scale_y_continuous(trans=scales::log2_trans(),breaks=c(1,3,5,10,100,1000))

p


