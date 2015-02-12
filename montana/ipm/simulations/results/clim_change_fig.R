##################################################
#### Script for figure of climate change results
#### Transfer to manuscript.rmd after complete
####
#### Andrew Tredennick
#### 2-11-2015
####

library(ggplot2)


files <- list.files()
cover_files <- files[grep("cover", files)]
spp_id <- substr(cover_files, 1, 4)
clim_id <- rep(c("observed", "pptChange", "tempChange", "temppptChange"), 4)
num_files <- length(cover_files)
all_sims <- data.frame(time=NA, cover=NA, species=NA, sim=NA)
for(i in 1:num_files){
  tmp <- read.csv(cover_files[i])
  tmp$species <- rep(spp_id[i], nrow(tmp))
  tmp$sim <- rep(clim_id[i], nrow(tmp))
  all_sims <- rbind(all_sims, tmp)
}
all_sims <- all_sims[2:nrow(all_sims),]

ggplot(all_sims)+
  geom_boxplot(aes(x=species, y=cover*100, fill=sim))+
  coord_cartesian(ylim = c(0,100))
