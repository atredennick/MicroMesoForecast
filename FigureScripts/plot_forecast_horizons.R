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
out_ipm <- list() # empty storage for data frames
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2ipms,tmp_file))
  tmp_sims$error <- with(tmp_sims, pred.cover.t1-obs.cover.t1)
  tmp_agg <- ddply(tmp_sims, .(horizon,species), summarise,
                   num_preds = length(error),
                   mean_error = mean(abs(error)),
                   sd_error = sd(abs(error)),
                   rmse = sqrt(mean(error^2)),
                   rho = cor(pred.cover.t1,obs.cover.t1))
  out_ipm <- rbind(out_ipm, tmp_agg)
}


##  QBM results
out_qbm <- list()
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2qbms,tmp_file))
  tmp_agg1 <- ddply(tmp_sims, .(horizon,species,quad,year.start,year.pred), summarise,
                    mean_pred = mean(pred.cover.t1),
                    mean_obs = mean(obs.cover.t1))
  tmp_agg1$error <- with(tmp_agg1, mean_pred - mean_obs)
  tmp_agg <- ddply(tmp_agg1, .(horizon,species), summarise,
                   num_preds = length(error),
                   mean_error = mean(abs(error)),
                   sd_error = sd(abs(error)),
                   rmse = sqrt(mean(error^2)),
                   rho = cor(mean_pred,mean_obs))
  out_qbm <- rbind(out_qbm, tmp_agg)
}

##  Combine the results
out_ipm$model <- "AIPM"
out_qbm$model <- "BQBM"
out_all <- rbind(out_ipm, out_qbm)



####
####  MAKE THE PLOT
####
out_plot <- subset(out_all, num_preds>49)
mycols <- c("#969C43", "#B46FA1")
mycols <- c("grey50", "black")
ggplot(data=out_plot, aes(x=horizon, y=rho))+
#   geom_ribbon(aes(x=yearsbefore, ymax=(mean_error+sd_error)*100, 
#                   ymin=(mean_error-sd_error)*100), alpha=0.3, color=NA)+
  geom_line(aes(linetype=model))+
  geom_point(aes(shape=model),size=2)+
  facet_grid(species~., scales="free")+
  xlab("Forecast Horizon (years)")+
  ylab(expression(paste("Forecast Accuracy (", rho, ")")))+
  # ylab("Mean Absolute Error (% Cover)")+
  scale_x_continuous(breaks=seq(1,13,by=2))+
  scale_shape_discrete(name="", labels=c("IPM", "QBM"))+
  scale_linetype_discrete(name="", labels=c("IPM", "QBM"))+
  theme_few()+
  theme(legend.position=c(0.8,0.96),
        legend.background=element_rect(NA))
ggsave(paste0(path2figs,"forecast_horizon.png"), width = 3, height = 6, dpi=120)
