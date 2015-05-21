###################################################
## Script to quick look at 1-step-ahead IPM output
##

spp_list <- c("BOGR", "HECO", "PASM", "POSE")
all_d <- data.frame(year=NA, lagcov=NA, obscov=NA, predcov=NA, 
                    group=NA, rep=NA, species=NA, quad=NA)
for(i in 1:length(spp_list)){
  tmp_d <- readRDS(paste("onestep_", spp_list[i],".RDS",sep=""))
#   tmp_d$species <- rep(spp_list[i],nrow(tmp_d))
  all_d <- rbind(all_d, tmp_d)
}
all_d <- all_d[2:nrow(all_d),]
all_d$resid <- with(all_d, (predcov*100)-(obscov*100))

saveRDS(all_d,"qbm_one-step_forecasts_combined.rds")

library(ggplot2)
ggplot(all_d, aes(x=obscov, y=predcov))+
  geom_point(alpha=0.2, shape=1)+
  geom_abline(aes(intercept=0, slope=1), color="red")+
  facet_wrap("species", scales="free")


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
               pearson_rho = cor(predcov, obscov),
               mean_cover = mean(obscov*100))
stats



