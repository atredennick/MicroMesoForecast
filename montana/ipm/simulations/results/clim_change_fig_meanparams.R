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
cover_files <- cover_files[grep("mean", cover_files)]
spp_id <- substr(cover_files, 1, 4)
clim_id <- rep(c("observed", "pptChange", "tempChange", "temppptChange"), 4)
num_files <- length(cover_files)
all_sims <- data.frame(time=NA, cover=NA, species=NA, sim=NA)
for(i in 1:num_files){
  tmp <- read.csv(cover_files[i])
#   print(nrow(tmp))
  tmp$species <- rep(spp_id[i], nrow(tmp))
  tmp$sim <- rep(clim_id[i], nrow(tmp))
  all_sims <- rbind(all_sims, tmp)
}
all_sims <- all_sims[2:nrow(all_sims),]

ggplot(all_sims)+
  geom_boxplot(aes(x=species, y=log(cover), fill=sim))
#   coord_cartesian(ylim = c(0,100))
tmp_sims <- subset(all_sims, cover < 1)
ggplot(tmp_sims, aes(x=time, y=cover*100, color=sim))+
  geom_line()+
  facet_wrap("species", scales="free")

simcols <- dcast(all_sims, species+time~sim, value.var = "cover")
obstmp <- simcols$observed
pertstmp <- simcols[,c("pptChange", "tempChange", "temppptChange")]
difffxn <- function(X){(log(X)-log(obstmp))}
tmpout <- as.data.frame(apply(pertstmp, MARGIN = 2, difffxn))
tmpout$species <- simcols$species
tmpm <- melt(tmpout, id.vars = "species")
all_means <- ddply(tmpm, .(species, variable), summarise,
                   avg_cover = mean(value),
                   med_cover = median(value),
                   sd_cover = sd(value),
                   up_cover = quantile(value, probs=0.875),
                   lo_cover = quantile(value, probs=0.125))
dgd = position_dodge(width = 0.60)
ggplot(all_means)+
  geom_point(aes(x=species, y=med_cover, shape=variable), 
             size=5, position=dgd)+
  geom_errorbar(aes(x=species, ymax=up_cover, ymin=lo_cover, group=variable),
                position=dgd, width=0.25)

saveRDS(all_means, "ipm_climatesims_logdiffs_mean.RDS")


# 
# tmpid <- which(all_means$sim=="observed")
# diffs <- list()
# for(i in 1:length(tmpid)){
#   obs <- all_means[tmpid[i],"avg_cover"]
#   diffs[[i]] <- (all_means[(tmpid[i]+1):(tmpid[i]+3), "avg_cover"] - rep(obs,3))/rep(obs,3)
# #   diffs[[i]] <- log(all_means[(tmpid[i]+1):(tmpid[i]+3), "avg_cover"] / rep(obs,3))
# }
# names(diffs) <- unique(all_sims$species)
# diff_df <- melt(as.data.frame(diffs))
# diff_df$sim <- rep(c("pptChange", "tempChange", "temppptChange"),4)
# 
# myCols2 <- c("grey45", "#277BA8", "#7ABBBD", "#AED77A")
# ggplot(diff_df, aes(x=variable, y=value, fill=sim))+
#   geom_bar(stat="identity", position="dodge", color="white")+
#   geom_hline(aes(yintercept=0))+
#   scale_fill_manual(values = myCols2[2:4])+
#   xlab("Species")+
#   ylab("Cover change (%)")+
#   theme_bw()
# 
