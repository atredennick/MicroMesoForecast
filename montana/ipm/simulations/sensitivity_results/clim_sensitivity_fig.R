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
clim_id <- rep(c("allClim", "Ppt", "Temp",
                 "Ppt", "Temp",
                 "Ppt", "Temp"), 4)
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
all_clim_sims <- subset(all_sims, vital=="All")

species <- unique(all_sims$species)
vital_adds <- c("Growth", "Survival", "Recruitment") 
out_controls <- data.frame(time=NA, cover=NA, species=NA, sim=NA, vital=NA)
for(spp in species){
  tmp <- subset(all_clim_sims, species==spp)
  for(vital in vital_adds){
    tmp$vital <- vital
    out_controls <- rbind(out_controls, tmp)
  }
}
out_controls <- out_controls[2:nrow(out_controls),]

all_noclim_sims <- subset(all_sims, vital!="All")
final <- rbind(out_controls, all_noclim_sims)

equilibrium_cover <- ddply(final, .(species, sim, vital), summarise,
                           eq_cover = median(cover))

mean_cover <- ddply(subset(equilibrium_cover, sim=="allClim"), .(species), summarise,
                    value = mean(eq_cover))

myCols2 <- c("#277BA8", "#7ABBBD", "#AED77A")
ggplot(equilibrium_cover, aes(x=sim, y=eq_cover*100, 
                              linetype=vital, group=vital, shape=vital))+
  geom_hline(data=mean_cover, aes(yintercept=value*100), linetype=2, color="grey45")+
  geom_line()+
  geom_point(size=4)+
  facet_grid(species~., scales = "free")+
  scale_linetype_manual(values=c(1,2,3), name="Vital Rate")+
  scale_shape_manual(values=c(19,17,15), name="Vital Rate")+
  ylab("Equilibrium Cover (%)")+
  xlab("Simulation")+
  theme_bw()
