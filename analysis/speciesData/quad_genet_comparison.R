rm(list=ls())
library(dplyr)
spp <- "BOGR"

quad_dat <- subset(read.csv("quadAllCover.csv"), Species=="BOGR")
gen_dat <- read.csv("BOGR/edited/survD.csv")
gen_based_cover <- ddply(gen_dat, .(quad,year), summarise,
                         gencover = sum(area)/10000)
combo <- merge(quad_dat, gen_based_cover, all=T)
combo[which(combo$propCover==0),"propCover"] <- NA
combo <- combo[which(combo$year!=45),]
combo$diff <- round(combo$propCover,4)-round(combo$gencover,4)
mean(combo[complete.cases(combo),"propCover"])
mean(combo[complete.cases(combo),"gencover"])
