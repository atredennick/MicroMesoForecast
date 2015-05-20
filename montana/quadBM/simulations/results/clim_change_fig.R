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
clim_id <- rep(c("noClimChange", "pptChange", "tempChange", "temppptChange"), 4)
num_files <- length(cover_files)
all_sims <- data.frame(time=NA, cover=NA, species=NA, sim=NA)
for(i in 1:num_files){
  tmp <- as.data.frame(readRDS(cover_files[i]))
  colnames(tmp) <- "cover"
  tmp$time <- seq(1:nrow(tmp))
  tmp$species <- rep(spp_id[i], nrow(tmp))
  tmp$sim <- rep(clim_id[i], nrow(tmp))
  tmp <- tmp[,c("time", "cover", "species", "sim")]
  all_sims <- rbind(all_sims, tmp)
}
all_sims <- all_sims[2:nrow(all_sims),]

ggplot(all_sims)+
  geom_boxplot(aes(x=species, y=cover*100, fill=sim))+
  coord_cartesian(ylim = c(0,100))

all_means <- ddply(all_sims, .(species, sim), summarise,
                   avg_cover = median(cover))
tmpid <- which(all_means$sim=="noClimChange")
diffs <- list()
for(i in 1:length(tmpid)){
  obs <- all_means[tmpid[i],"avg_cover"]
  diffs[[i]] <- (all_means[(tmpid[i]+1):(tmpid[i]+3), "avg_cover"] - rep(obs,3))/rep(obs,3)
}
names(diffs) <- unique(all_sims$species)
diff_df <- melt(as.data.frame(diffs))
diff_df$sim <- rep(c("pptChange", "tempChange", "temppptChange"),4)
saveRDS(diff_df, "qbm_climatesims_percdiffs.rds")
myCols2 <- c("grey45", "#277BA8", "#7ABBBD", "#AED77A")
ggplot(diff_df, aes(x=variable, y=value, fill=sim))+
  geom_bar(stat="identity", position="dodge", color="white")+
  geom_hline(aes(yintercept=0))+
  scale_fill_manual(values = myCols2[2:4])+
  xlab("Species")+
  ylab("Cover change (%)")+
  theme_bw()

