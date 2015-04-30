##  Script to look at results from leave-one-year-out growth fits

rm(list=ls(all=TRUE))

library(plyr)
library(reshape2)
library(ggplot2)

all_files <- list.files()
rds_files <- all_files[grep("*.rds", all_files)]

test <- melt(readRDS(rds_files[2]))
numvars <- length(unique(test$Var2))
ggplot(test, aes(x=Var1, y=value))+
  geom_line()+
  facet_wrap("Var2", scales = "free_y")