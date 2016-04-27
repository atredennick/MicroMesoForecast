##  R Script to Look at Sensitivity of Individual Vital Rates to 
##  Perturbed Climate Covariates as Applied to an Average-Sized Plant

##  Author: Andrew Tredennick
##  Last update: 4-27-2016

# Clear the workspace
rm(list=ls())



####
####  PRELIMINARIES
####
path2genets <- "../analysis/data_processing/speciesData/"
path2ipms <- "../analysis/ipm/vitalRateRegs/growth/"
species_names <- c("BOGR", "HECO", "PASM", "POSE")



####
####  DEFINE SIMPLE GROWTH FUNCTION FOR A SINGLE PLANT
####
grow_plant <- function(x, intercept, slope, climate_effects, weather_vector) {
  mu <- intercept + slope*x + sum(climate_effects*weather_vector)
  return(mu)
}



####
####  READ IN MCMC AND BREAK APART
####
ispp <- 1
do_species <- species_names[ispp]
posteriors <- readRDS(paste0(path2ipms, "growth_stanmcmc_", do_species))





