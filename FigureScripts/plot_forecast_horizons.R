##  R Script to plot forecast horizons for both model types and all species.
##
##  Author: Andrew Tredennick
##  Last update: 4-20-2016

# Clear the workspace
rm(list=ls())


####
####  PRELIMINARIES
####
path2ipms <- "../analysis/ipm/simulations/results/forecast_horizon/"
path2qbms <- "../analysis/quadBM/simulations/results/forecast_horizon/"
path2figs <- "../manuscript/components/"
species_list <- c("BOGR", "HECO", "PASM", "POSE")



####
####  LOAD LIBRARIES
####
library(ggplot2)  # plotting
library(ggthemes) # plotting themes
library(reshape2) # data wrangling
library(plyr)     # data wrangling/summarizing



####
####  READ IN SIMULATION RESULTS AND COLLATE/FORMAT
####
# Saved simulation results are the final year (t) of cover predicted by the 
# model from some starting year t-x, where x ranges from 1 to 12 years. Each
# forecast simulation was run 50 times with different parameters so that we can
# plot uncertainty around the forecast horizon trend.

##  IPM results
out_ipm <- data.frame(startyear=NA, mean_error=NA, rmse=NA, species=NA)
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2ipms,tmp_file))
  tmp_sims$error <- with(tmp_sims, finalyear_cover-obs.cover.t1)
  tmp_agg <- ddply(tmp_sims, .(startyear), summarise,
                   mean_error = mean(error),
                   rmse = sqrt(mean(error^2)))
  tmp_agg$species <- species_list[ispp]
  out_ipm <- rbind(out_ipm, tmp_agg)
}
out_ipm <- out_ipm[2:nrow(out_ipm),] # remove first NA-filled row

##  QBM results
out_qbm <- data.frame(yearstart=NA, mean_error=NA, rmse=NA, species=NA)
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2qbms,tmp_file))
  tmp_sims <- tmp_sims[2:nrow(tmp_sims),]
  tmp_sims$error <- with(tmp_sims, finalcover-obs_finalcover)
  tmp_agg <- ddply(tmp_sims, .(yearstart), summarise,
                   mean_error = mean(error),
                   rmse = sqrt(mean(error^2)))
  tmp_agg$species <- species_list[ispp]
  out_qbm <- rbind(out_qbm, tmp_agg)
}
out_qbm <- out_qbm[2:nrow(out_qbm),] # remove first NA-filled row

##  Combine the results
out_ipm$model <- "IPM"
out_ipm$yearsbefore <- with(out_ipm, max(startyear)-startyear+1)
out_qbm$model <- "QBM"
out_qbm$yearsbefore <- with(out_qbm, max(yearstart)-yearstart+1)
names(out_qbm) <- names(out_ipm)
out_all <- rbind(out_ipm, out_qbm)



####
####  MAKE THE PLOT
####
ggplot(data=out_all, aes(x=yearsbefore, y=abs(mean_error)*100, linetype=model, shape=model))+
  geom_line(color="grey50")+
  geom_point(size=3)+
  facet_grid(species~., scales="free")+
  xlab("Years Before Forecast")+
  # ylab(expression(atop("Forecast Skill", atop("(1 - Root Mean Square Error)", ""))))+
  ylab("Mean Absolute Error (% Cover)")+
  scale_x_continuous(breaks=c(1:13))+
  scale_shape_discrete(name="")+
  scale_linetype_discrete(name="")+
  theme_few()+
  theme(legend.position=c(0.85,0.97),
        legend.background=element_rect(NA))