rm(list = ls())
library(pophelper)

afiles <- list.files("./input_for_pop_k10/", pattern = "M.+Q$", full.names = T)
# create a qlist
alist <- readQ(afiles)
col <- c(
  "#002cf9", "#ff2920", "#00fb3a", "#ff942d",
  "#9b1e8f", "#ffff44", "#FDF060", "#A945FF",
  "#BFF217", "#60D5FD", "#CC1577", "#F2B950",
  "#7FB21D", "#EC496F", "#326397", "#B26314",
  "#027368", "#A4A4A4", "#610B5E", "#a6cee3",
  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c"
)



g <- read.table("./Group_fastSt.txt",
                head = T, stringsAsFactors = F
)

plotQ(alist[c(1:12)],
      imgoutput = "join",
      clustercol = col,
      grplab = g,
      panelratio = c(1, 1),
      grplabsize = 0.8, linesize = 0.1, pointsize = 0.5,
      grplabpos = 0.9, grplabangle = 90,
      ordergrp = TRUE,
      subsetgrp = c(
        "DOM_CID",
        "DOM_DES",
        "SYL_DNK",
        "SYL_FRA_WEST",
        "SYL_FRA_EAST",
        "SYL_DEU",
        "SYL_AUT",
        "SYL_ROU",
        "SYL_MKD",
        "SYL_ITA",
        "ORI_ARM",
        "SIE_KAZ",
        "SIE_CHN",
        "BAC_CHN"
      ),
      divsize = 0.1,
      indlabheight = 0,
      width = 14,
      height = 1,
      panelspacer = 0.01,
      sppos = "left",
      imgtype = "pdf",
      linepos = 0.95,
      exportpath = getwd(),
      splabsize = 4,
      splab = c(
        "K=2 Major",
        "K=3 Major",
        "K=4 Major",
        "K=5 Major",
        "K=6 Major",
        "K=6 Minor",
        "K=7 Major",
        "K=8 Major",
        "K=9 Major",
        "K=9 Minor",
        "K=10 Major",
        "K=10 Minor"
      )
)

