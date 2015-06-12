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
num_files <- length(cover_files)
all_sims <- data.frame(time=NA, cover=NA, species=NA, climsim=NA, vital=NA)
for(i in 1:num_files){
  tmp <- read.csv(cover_files[i])
  tmp$species <- rep(spp_id[i], nrow(tmp))
  all_sims <- rbind(all_sims, tmp)
}
all_sims <- all_sims[2:nrow(all_sims),]

# Rename vital rates with codes
vital_names <- unique(all_sims$vital)
vital_codes <- c("A", "G", "GR", "GS", "R", "S" , "SR") 
vitals <- data.frame(names=vital_names,code=vital_codes)
for(vnow in vital_names){
  ids <- which(all_sims[,"vital"]==vnow)
  idcode <- which(vitals[,"names"]==vnow)
  all_sims[ids,"vital"] <- as.character(vitals[idcode,"code"])
}

# Rep the 'all' simulation for each vital rate for plotting
all_clim_sims <- subset(all_sims, vital=="A")
species <- unique(all_sims$species)
out_controls <- data.frame(time=NA, cover=NA, species=NA, climsim=NA, vital=NA)
for(spp in species){
  tmp <- subset(all_clim_sims, species==spp)
  for(vital in vital_codes){
    tmp$vital <- vital
    out_controls <- rbind(out_controls, tmp)
  }
}
out_controls <- out_controls[2:nrow(out_controls),]

# Combine the datasets
all_noclim_sims <- subset(all_sims, vital!="all")
final <- rbind(out_controls, all_noclim_sims)

# Calculate statistics for plotting
equilibrium_cover <- ddply(final, .(species, climsim, vital), summarise,
                           eq_cover = median(cover))
mean_cover <- ddply(subset(equilibrium_cover, climsim=="all"), .(species), summarise,
                    value = mean(eq_cover))

# myCols2 <- c("#277BA8", "#7ABBBD", "#AED77A")
ggplot(equilibrium_cover, aes(x=climsim, y=eq_cover*100, 
                              color=vital, group=vital))+
  geom_hline(data=mean_cover, aes(yintercept=value*100), linetype=2, color="grey45")+
  geom_line()+
  geom_point(size=4)+
  facet_grid(species~., scales = "free")+
#   scale_linetype_manual(values=c(1,2,3), name="Vital Rate")+
#   scale_shape_manual(values=c(19,17,15), name="Vital Rate")+
  ylab("Equilibrium Cover (%)")+
  xlab("Simulation")+
  theme_bw()
