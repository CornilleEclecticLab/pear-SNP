

# This script refers to the following files:
# https://github.com/Jimi92/Population-genomics-Pyrenophora-teres/blob/main/selection/get_selection_figure.Rmd

rm(list=ls())
library(ggplot2)

setwd('/shared/home/ynie/work/pear/run_OmegaPlus/bin')



read_omega_pops <- function(f,pop) {
  df <- read.table(f, skip = 2, header = FALSE)
  colnames(df) <- c("Position", "Omega","Left_boundary", "Right_boundary","Validity")
  mm = median(df$Omega)
  df["standar_omega"] = df["Omega"]/mm
  df["Population"] <- pop
  return(df)
}

c1 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min1000", "EW_Grid10K_min1K")
c2 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min2000", "EW_Grid10K_min2K")
c3 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min2500", "EW_Grid10K_min2.5K")
c4 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min4000", "EW_Grid10K_min4K")
c5 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min5000", "EW_Grid10K_min5K")


c1 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid1000.min1000", "EW_Grid1K_min2.5K")
c2 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid2000.min2000", "EW_Grid2K_min2.5K")
c3 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid4000.min2500", "EW_Grid4K_min2.5K")
c4 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid6000.min4000", "EW_Grid6K_min2.5K")
c5 <- read_omega_pops("./test/OmegaPlus_Report.s02.run_OmegaPlus.EuropWild.GWHBAOS00000076.grid10000.min5000", "EW_Grid10K_min2.5K")




full_df <- rbind(c2,c1)
full_df <- rbind(full_df,c3)
full_df <- rbind(full_df,c4)
full_df <- rbind(full_df,c5)
full_df <- rbind(full_df,c5)

library(ggplot2)

ggplot(full_df,aes(Position,Omega, color = Population )) + 
  geom_line() + 
  # geom_hline(yintercept = 36, colour = "red") + 
  ggtitle("Chromosome 1") + 
  xlim(0,2900000) + 
  theme_minimal() 

ggplot(full_df,aes(Position,standar_om, color = Population)) + 
  geom_line (a=0.3, size = 1) + ggtitle("                                               Chromosome 6") + 
  xlim(2500000,2900000) + 
  ylim(0,20)  + 
  theme_minimal()



##################################################################################################
pdf("ND_chr6.pdf")
ggplot(c1,aes(Position,standar_om, color = standar_om < 2.5)) + 
  geom_point(alpha=0.3) + 
  ggtitle("Chromosome 6, North Dakota") + 
  xlim(2860000,2900000) + 
  ylim(0,20)  + 
  geom_vline(xintercept = 2874993) + 
  scale_colour_manual(name = 'Significant ', values = setNames(c('grey','aquamarine4'),c(T, F))) + 
  geom_vline(xintercept = 2879103)+ 
  ylab("Omega") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
dev.off

pdf('Moroc_chr6.pdf')
ggplot(c3,aes(Position,standar_om, color = standar_om < 2.5)) + geom_point(alpha=0.3) + ggtitle("Chromosome 6, Morocco") + xlim(2860000,2900000) + ylim(0,20) +
  scale_colour_manual(name = 'Significant ', values = setNames(c('grey','aquamarine4'),c(T, F))) +
  geom_vline(xintercept = 2874993) + 
  geom_vline(xintercept = 2879103)+ ylab("Omega") + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
dev.off

pdf('France_chr6.pdf')
ggplot(c4,aes(Position,standar_om, color = standar_om < 5.2)) + geom_point(alpha=0.3) + ggtitle("Chromosome 6, France") + xlim(2860000,2900000) + ylim(0,20) + geom_vline(xintercept = 2874993) + 
  scale_colour_manual(name = 'Significant ', values = setNames(c('grey','aquamarine4'),c(T, F))) +
  geom_vline(xintercept = 2879103)+ ylab("Omega") + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
dev.off

pdf('Azer_chr6.pdf')
ggplot(c5,aes(Position,standar_om, color = standar_om < 2.5)) + geom_point(alpha=0.3) + 
  scale_colour_manual(name = 'Significant ', values = setNames(c('grey','aquamarine4'),c(T, F))) +
  ggtitle("Chromosome 6, Azerbaijan") + xlim(2860000,2900000) + ylim(0,20) + geom_vline(xintercept = 2874993) + geom_vline(xintercept = 2879103)+ ylab("Omega") + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))
dev.off
