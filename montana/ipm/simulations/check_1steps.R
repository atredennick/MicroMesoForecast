###################################################
## Script to quick look at 1-step-ahead IPM output
##

file <- "_sim_cover_1step_ahead_new.csv"
spp_list <- c("BOGR", "HECO", "PASM", "POSE")
all_d <- data.frame(quad=NA,t1=NA,t0=NA,rep=NA,cover.t0=NA,cover.t1=NA,obs.cover.t1=NA,species=NA)
for(i in 1:length(spp_list)){
  tmp_d <- read.csv(paste(spp_list[i],file,sep=""))
  tmp_d$species <- rep(spp_list[i],nrow(tmp_d))
  all_d <- rbind(all_d, tmp_d)
}
all_d <- all_d[2:nrow(all_d),]
all_d$resid <- with(all_d, (obs.cover.t1*100)-(cover.t1*100))

library(ggplot2)
ggplot(all_d,aes(x=t1,y=resid,group=t1))+
  geom_boxplot()+
  facet_wrap("species", scale="free")




