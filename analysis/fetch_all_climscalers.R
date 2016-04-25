##  Script to collate climate covariate mean and sds for scaling
##  Calculates scalers by species, vital rate, and left-out year

# Clear the workspace
rm(list=ls())

species <- c("BOGR", "HECO", "PASM", "POSE")
clim_covariates <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", 
                     "TmeanSpr2", "ppt1TmeanSpr1", "ppt2TmeanSpr2")

####
####  PERCENT COVER COVARIATES
####
##  Read in data
#bring in data
cover_dat <- readRDS("./processed_data/cover_with_weather.RDS")
years <- c(1,unique(cover_dat$year))
full_df <- list()
for(do_species in species){
  for(rm_year in 1:length(years)){
    growD <- subset(cover_dat, species==do_species & year!=years[rm_year])
    clim_covs <- growD[,clim_covariates]
    # Get scalers for climate covariates from training data
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("qbm", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df, "qbm_all_clim_scalers.RDS")



####
####  GROWTH PROCESS
####
grow_dat <- readRDS("./processed_data/grow_with_weather.RDS")
years <- c(1,unique(grow_dat$year))
full_df <- list()
for(do_species in species){
  for(rm_year in 1:length(years)){
    growD <- subset(grow_dat, species==do_species & year!=years[rm_year])
    clim_vars_all <- c(clim_covariates)
    clim_covs <- growD[,clim_vars_all]
    # Get scalers for climate covariates from training data
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("growth", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df, "growth_all_clim_scalers.RDS")



####
####  SURVIVAL PROCESS
####
surv_dat <- readRDS("./processed_data/surv_with_weather.RDS")
years <- c(1,unique(surv_dat$year))
full_df <- list()
for(do_species in species){
  for(rm_year in 1:length(years)){
    survD <- subset(surv_dat, species==do_species & year!=years[rm_year])
    clim_covs <- survD[,clim_covariates]
    # Get scalers for climate covariates from training data
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("survival", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df, "survival_all_clim_scalers.RDS")



####
####  RECRUITMENT PROCESS
####
rec_dat <- readRDS("./processed_data/rec_with_weather.RDS")
years <- c(1,unique(rec_dat$year))
full_df <- list()
for(do_species in species){
  for(rm_year in 1:length(years)){
    recD <- subset(rec_dat, species==paste("R.",do_species,sep="") & year!=years[rm_year])
    clim_covs <- recD[,clim_covariates]
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("recruitment", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df, "recruitment_all_clim_scalers.RDS")
