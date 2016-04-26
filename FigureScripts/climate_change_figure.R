##  R Script for Plotting Results from Climate Change Scenarios
##
##  Author: Andrew Tredennick
##  Last update: 4-26-2016

# Clear the workspace
rm(list=ls())



####
####  PRELIMINARIES
####
path2qbms <- "../analysis/quadBM/simulations/results/climchange_varyparams/"
path2ipms <- "../analysis/ipm/simulations/results/climchange_varyparams/"
path2figs <- "../manuscript/components/"



####
####  LOAD LIBRARIES
####
library("plyr") # data summarising
library("reshape2") # data wrangling
library("ggplot2") # plotting
library("ggthemes") # plotting themes



####
####  LOAD AND CONFIGURE QBM RESULTS
####
qbm_files <- list.files(path2qbms)
qbm_allsims <- list()
for(ifile in 1:length(qbm_files)){
  tmpfile <- qbm_files[ifile]
  tmp <- read.csv(paste0(path2qbms, tmpfile))
  tmp_sim <- strsplit(strsplit(tmpfile, "_")[[1]][4], "[.]")[[1]][1]
  tmp$climsim <- tmp_sim
  qbm_allsims <- rbind(qbm_allsims, tmp)
}
# Test for proper number of climate change simulations
if(length(unique(qbm_allsims$climsim))!=4) { stop("wrong number of climate perturbations") }

## Calculate average and +/- proportional differences for each clim perturbation
# First, add in no perturbation climate results as new column, by species
species_names <- unique(qbm_allsims$species)
qbm_diffs <- list()
for(ispp in 1:length(species_names)){
  tmp <- subset(qbm_allsims, species==species_names[ispp])
  obsclim_cover <- subset(tmp, climsim=="NAChange")["cover"]
  tmp2 <- subset(tmp, climsim!="NAChange")
  tmp2$obsclim_cover <- rep(obsclim_cover$cover, times = length(unique(tmp2$climsim)))
  qbm_diffs <- rbind(qbm_diffs, tmp2)
}
qbm_diffs$propdiff <- with(qbm_diffs, log(cover*100)-log(obsclim_cover*100))

# Calculate mean proportional difference by species and climate perturbation
# and associated 90% quantiles
qbm_diffs_aggregated <- ddply(qbm_diffs, .(climsim, species), summarise,
                              mean_prop_diff = mean(propdiff),
                              upper_diff = quantile(propdiff, probs = 0.95),
                              lower_diff = quantile(propdiff, probs = 0.05))



####
####  LOAD AND CONFIGURE IPM RESULTS
####
ipm_files <- list.files(path2ipms)
ipm_allsims <- list()
for(ifile in 1:length(ipm_files)){
  tmpfile <- ipm_files[ifile]
  tmp <- read.csv(paste0(path2ipms, tmpfile))
  tmp_sim <- strsplit(strsplit(tmpfile, "_")[[1]][4], "[.]")[[1]][1]
  tmp$climsim <- tmp_sim
  ipm_allsims <- rbind(ipm_allsims, tmp)
}
# Test for proper number of climate change simulations
if(length(unique(ipm_allsims$climsim))!=4) { stop("wrong number of climate perturbations") }

## Calculate average and +/- proportional differences for each clim perturbation
# First, add in no perturbation climate results as new column, by species
species_names <- unique(ipm_allsims$species)
ipm_diffs <- list()
for(ispp in 1:length(species_names)){
  tmp <- subset(ipm_allsims, species==species_names[ispp])
  obsclim_cover <- subset(tmp, climsim=="NAChange")["cover"]
  tmp2 <- subset(tmp, climsim!="NAChange")
  tmp2$obsclim_cover <- rep(obsclim_cover$cover, times = length(unique(tmp2$climsim)))
  ipm_diffs <- rbind(ipm_diffs, tmp2)
}
ipm_diffs$propdiff <- with(ipm_diffs, log(cover*100)-log(obsclim_cover*100))

# Calculate mean proportional difference by species and climate perturbation
# and associated 90% quantiles
ipm_diffs_aggregated <- ddply(ipm_diffs, .(climsim, species), summarise,
                              mean_prop_diff = mean(propdiff),
                              upper_diff = quantile(propdiff, probs = 0.95),
                              lower_diff = quantile(propdiff, probs = 0.05))



####
####  COMBINE RESULTS AND PLOT
####
ipm_diffs_aggregated$model <- "ipm"
qbm_diffs_aggregated$model <- "qbm"
diffs_aggregated <- rbind(ipm_diffs_aggregated, qbm_diffs_aggregated)

# Change some names for the plot
diffs_aggregated[which(diffs_aggregated$climsim=="pptChange"), "climsim"] <- "Appt"
diffs_aggregated[which(diffs_aggregated$climsim=="TmeanChange"), "climsim"] <- "Btmean"
diffs_aggregated[which(diffs_aggregated$climsim=="pptTmeanChange"), "climsim"] <- "Cboth"

dgd <- position_dodge(width = 0.6)
ggplot(diffs_aggregated, aes(x=climsim, y=mean_prop_diff))+
  geom_hline(aes(yintercept=0), color="grey", linetype=2)+
  geom_errorbar(aes(ymin=lower_diff, ymax=upper_diff, group=model), width=0.2, position=dgd)+
  geom_point(aes(shape=model), position=dgd, size=4)+
  facet_wrap("species", ncol=1, scales="free_y")+
  scale_x_discrete(labels=c("+ppt", "+temp", "+ppt&temp"))+
  # scale_linetype_discrete(name="Model", labels=c("IPM", "QBM"))+
  scale_shape_discrete(name="Model", labels=c("IPM", "QBM"))+
  ylab("Proportional Change in Cover")+
  xlab("Climate Simulation")+
  theme_few()
ggsave(paste0(path2figs, "climate_change_results.png"), width = 4.5, height = 7, units = "in", dpi = 100)


