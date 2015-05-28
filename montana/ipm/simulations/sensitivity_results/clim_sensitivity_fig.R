##################################################
#### Script for figure of climate change results
#### Transfer to manuscript.rmd after complete
####
#### Andrew Tredennick
#### 2-11-2015
####

library(ggplot2)
library(plyr)
library(reshape2)


files <- list.files()
cover_files <- files[grep("cover", files)]
spp_id <- substr(cover_files, 1, 4)
clim_id <- rep(c("allClim", "noPpt", "noTemp",
                 "noPpt", "noTemp",
                 "noPpt", "noTemp"), 4)
vital_id <- rep(c("All", "Growth", "Growth",
                  "Survival", "Survival",
                  "Recruitment", "Recruitment"), 4)
num_files <- length(cover_files)
all_sims <- data.frame(time=NA, cover=NA, species=NA, sim=NA, vital=NA)
for(i in 1:num_files){
  tmp <- read.csv(cover_files[i])
  tmp$species <- rep(spp_id[i], nrow(tmp))
  tmp$sim <- rep(clim_id[i], nrow(tmp))
  tmp$vital <- rep(vital_id[i], nrow(tmp))
  all_sims <- rbind(all_sims, tmp)
}
all_sims <- all_sims[2:nrow(all_sims),]

ggplot(all_sims)+
  geom_boxplot(aes(x=species, y=cover, fill=sim))
#   coord_cartesian(ylim = c(0,100))

all_means <- ddply(all_sims, .(species, sim), summarise,
                   avg_cover = median(cover))

myCols2 <- c("grey45", "#277BA8", "#7ABBBD", "#AED77A")
ggplot(diff_df, aes(x=sim, y=value, fill=sim))+
  geom_bar(stat="identity", position="dodge", color="white")+
  geom_hline(aes(yintercept=0))+
  scale_fill_manual(values = myCols2[2:4])+
  xlab("Species")+
  ylab("Cover change (%)")+
  theme_bw()

