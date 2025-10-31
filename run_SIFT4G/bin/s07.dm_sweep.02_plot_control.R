# By Xilong CHEN
# Create date: 2023-09-21
# Contact: chen_xilong@outlook.com

# Set Working Directory to source file folder

rm(list = ls())

library("ggplot2")
library("cowplot")
library(ggpubr)


### read data
df_sweep <- read.table("../output/s03.dm_sweep/region.mutation_burden.Rdata.txt", header =T)

### use all of the dm
df_sweep_all <- df_sweep[df_sweep$type == "All", ]

#### plot five pops 

df_sweep_all$pop <- factor(df_sweep_all$pop, c("DomC", "DomD", "SiOr", "SylW", "SylE"))

df_sweep_all$region_type <- factor(df_sweep_all$region_type, c("control", "hard_sweep", "soft_sweep"))

df_sweep_all$de_per_100k_cds <- df_sweep_all$de_per_cds * 100000

p <- ggboxplot(df_sweep_all, x = "pop", y = "de_per_100k_cds",
               color = "region_type", palette = "npg"
) + theme_bw(base_size = 8)
# Use only p.format as label. Remove method name.
p + stat_compare_means(aes(group = region_type), label = "p.signif")

ggsave2("../output/s03.dm_sweep/pop_all.region.mutation_burden.svg", width=7, height=2.5)
ggsave2("../output/s03.dm_sweep/pop_all.region.mutation_burden.pdf", width=7, height=2.5)



### plot only DomC and DomD
df_sweep_all_dom <- df_sweep_all[df_sweep_all$pop == "DomC" | df_sweep_all$pop == "DomD", ]

df_sweep_all_dom$pop <- factor(df_sweep_all_dom$pop, c("DomC", "DomD"))

df_sweep_all_dom$region_type <- factor(df_sweep_all_dom$region_type, c("control", "hard_sweep", "soft_sweep"))

df_sweep_all_dom$de_per_100k_cds <- df_sweep_all_dom$de_per_cds * 100000

p <- ggboxplot(df_sweep_all_dom, x = "pop", y = "de_per_100k_cds",
               color = "region_type", palette = "npg"
) + theme_bw(base_size = 8)
# Use only p.format as label. Remove method name.
p + stat_compare_means(aes(group = region_type), label = "p.signif")
p + theme(legend.position='none')
ggsave2("../output/s03.dm_sweep/pop_Dom.region.mutation_burden.svg", width=2, height=2.5)
ggsave2("../output/s03.dm_sweep/pop_Dom.region.mutation_burden.pdf", width=2, height=2.5)

