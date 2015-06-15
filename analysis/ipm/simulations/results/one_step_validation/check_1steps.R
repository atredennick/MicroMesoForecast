###################################################
## Script to quick look at 1-step-ahead IPM output
##

file <- "_sim_cover_1step_ahead_year.csv"
spp_list <- c("BOGR", "HECO", "PASM", "POSE")
all_d <- data.frame(quad=NA,t1=NA,t0=NA,rep=NA,cover.t0=NA,cover.t1=NA,obs.cover.t1=NA,species=NA)
for(i in 1:length(spp_list)){
  tmp_d <- read.csv(paste(spp_list[i],file,sep=""))
  tmp_d$species <- rep(spp_list[i],nrow(tmp_d))
  all_d <- rbind(all_d, tmp_d)
}
all_d <- all_d[2:nrow(all_d),]
all_d$resid <- with(all_d, (cover.t1*100)-(obs.cover.t1*100))
all_d$year <- rep("ayear effect", nrow(all_d))

# file <- "_sim_cover_1step_ahead_noYear.csv"
# spp_list <- c("BOGR", "HECO", "PASM", "POSE")
# all_d2 <- data.frame(quad=NA,t1=NA,t0=NA,rep=NA,cover.t0=NA,cover.t1=NA,obs.cover.t1=NA,species=NA)
# for(i in 1:length(spp_list)){
#   tmp_d <- read.csv(paste(spp_list[i],file,sep=""))
#   tmp_d$species <- rep(spp_list[i],nrow(tmp_d))
#   all_d2 <- rbind(all_d2, tmp_d)
# }
# all_d2 <- all_d2[2:nrow(all_d2),]
# all_d2$resid <- with(all_d2, (obs.cover.t1*100)-(cover.t1*100))
# all_d2$year <- rep("no year effect", nrow(all_d2))

# all_d <- rbind(all_d, all_d2)

saveRDS(all_d,"ipm_loyo_forecasts_combined.rds")

library(ggplot2)
library(ggplot2)
ggplot(all_d, aes(x=obs.cover.t1, y=cover.t1))+
  geom_point()+
  geom_abline(aes(intercept=0, slope=1), color="red")+
  facet_wrap("species", scales="free")
# 
# myCols <- c("#237DA4", "#7ECAD0")
# ggplot(all_d,aes(x=as.character(t1),y=resid))+
#   geom_boxplot(outlier.size = 0)+
#   facet_wrap("species", scale="free") #+
#   scale_fill_manual(values=myCols, labels=c("Year effect", "No year effect"))+
#   theme_bw()

#Calculate average residual variation by species
library(plyr)
stats <- ddply(all_d, .(species), summarise,
               mean_abs_error = mean(abs(resid)),
               pearson_rho = cor(cover.t1, obs.cover.t1),
               mean_cover = mean(obs.cover.t1*100))
stats



