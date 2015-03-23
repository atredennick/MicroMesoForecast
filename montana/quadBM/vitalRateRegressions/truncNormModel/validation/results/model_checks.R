##  Quick script to look at model checks for validation models

library('reshape2'); library('ggplot2')

## Gelman plots
gelman_files <- list.files()[grep("growthGelman*", list.files())]
gelman_list <- list()
for(gelmans in 1:length(gelman_files)){
  gelman_list[[gelmans]] <- read.csv(gelman_files[gelmans])
}
gelman_data <- melt(gelman_list)
ggplot(subset(gelman_data, variable=="Point.est."))+
  geom_histogram(aes(x=value), color="grey")

ggplot(subset(gelman_data, variable=="Point.est."))+
  geom_histogram(aes(x=value))+
  facet_wrap("X")