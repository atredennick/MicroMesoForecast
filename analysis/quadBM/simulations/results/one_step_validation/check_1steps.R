###################################################
## Script to quick look at 1-step-ahead IPM output
##

spp_list <- c("BOGR", "HECO", "PASM", "POSE")
all_d <-  list()
for(i in 1:length(spp_list)){
  tmp_d <- readRDS(paste(spp_list[i],"_sim_cover_1step_ahead_year.RDS",sep=""))
#   tmp_d$species <- rep(spp_list[i],nrow(tmp_d))
  all_d <- rbind(all_d, tmp_d)
}

# all_d$resid <- with(all_d, (predcov*100)-(obscov*100))
# removes <- which(is.na(all_d$sim)==TRUE)
# all_d <- all_d[-removes, ]
saveRDS(all_d,"qbm_one-step_forecasts_combined.rds")
# 
# library(ggplot2)
# ggplot(all_d, aes(x=cover.t0, y=cover.t1))+
#   geom_point(alpha=0.2, shape=1)+
#   geom_abline(aes(intercept=0, slope=1), color="red")+
#   facet_wrap("species", scales="free")


# myCols <- c("#237DA4", "#7ECAD0")
# ggplot(all_d,aes(x=as.character(t1),y=resid))+
#   geom_boxplot(outlier.size = 0)+
#   facet_wrap("species", scale="free") #+
#   scale_fill_manual(values=myCols, labels=c("Year effect", "No year effect"))+
#   theme_bw()

#Calculate average residual variation by species
all_d$resid <- with(all_d, (pred_cover.t1*100)-(obs_cover.t1*100))
library(plyr)
stats <- ddply(all_d, .(species), summarise,
               mean_abs_error = mean(abs(resid)))
stats



