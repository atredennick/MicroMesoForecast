##  R Script to Run Climate Change Scenario Simulations Using the
##  Quad-Based Model of Percent Cover Change.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Last update: April 11, 2016
##
##  Climate change scenarios are %age increases in:
##    1. Precipitation
##    2. Temperature
##    3. Precipitation and Temperature
##  where the %age increase is set by the user, but we focus on 1%
##  increases in mean climate variables for the paper. 

##  We define two functions to conduct these simulations: a cover change function 
##  based on the truncated log-normal statistical model of change in species'
##  proportional cover (sourced), and a climate change function that alter the climate
##  time series and standardizes the altered time series using mean and 
##  standard deviation from the observed climate time series as repeated in
##  the percent cover dataframe (e.g., multiple records of year 1933 climate).



# ----------------------------- ### CODE STARTS HERE ### ---------------------- #
##  Clear the Workspace
rm(list=ls())

##  Set working directory programmatically
root <- ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos") # modify as needed
setwd(paste(root,"/MicroMesoForecast/analysis/quadBM/simulations",sep="")) # modify as needed 


####
####  LOAD LIBRARIES
####
library(reshape2) # data wrangling
library(plyr) # data wrangling
library(EnvStats) # for truncated log-normal functions



####
####  PRELIMINARIES
####
##  Use mean parameter values?
use_mean_params <- FALSE

tsims <- 2500 # number of iterations ("years") to run the model
burn.in <- 500 # number of "years" to throw out from transient phase
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)

##  Degree of climate change (can be a vector)
# clim_change <- c(0.01, 0.1, 0.2, 0.3) # climate change percentage
do_change <- 1
clim_change <- 0.01 # climate change percentage in decimal form
perc_change <- clim_change*100 # climate change percentage in %age form
clim_alt <- clim_change[do_change] # can be set programmatically if necessary
filetag <- paste(perc_change[do_change], ".csv", sep="")

##  Load climate scalers
qbm_clim_scalers <- readRDS("../../qbm_all_clim_scalers.RDS")

##  Directory paths
if(use_mean_params==TRUE){
  resultspath <- "./results/climchange_meanparams/"
}
if(use_mean_params==FALSE){
  resultspath <- "./results/climchange_varyparams/"
}

climvec <- readRDS("../../climate_year_sequence.rds") # random climate year sequence
path2mcmcs <- "../vitalRateRegressions/truncNormModel/"


####
####  LOAD CLIMATE DATA (UNSTANDARDIZED)
####
clim_data_obs <- read.csv("../../weather/Climate.csv")
# Subset and reorder to match regression param import
clim_data_obs <- clim_data_obs[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")]
clim_covariates <- c("pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2","ppt1TmeanSpr1","ppt2TmeanSpr2")


####
####  SOURCE FUNCTION FILES
####
source("fetch_params_fxn.R")
source("qbm_sim_fxn.R")



####
####  FUNCTION FOR PERTURBING CLIMATE FOR IPM SIMULATIONS
####  AND FOR RETREIVING SPECIES-SPECIFIC CLIMATE SCALERS 
####
perturb_climate <- function(clim_alt, clim_var, clim_data_obs, do_spp){
  doSpp <- spp_list[do_spp]
  
  if(is.na(clim_var)==TRUE){
    clim_data <- clim_data_obs
    clim_data$ppt1TmeanSpr1 <- with(clim_data, ppt1*TmeanSpr1)
    clim_data$ppt2TmeanSpr2 <- with(clim_data, ppt1*TmeanSpr2)
  }
  
  if(is.na(clim_var)==FALSE){
    # Perturb climate time series
    clim_data <- clim_data_obs # re-assigns the temporary clim_data object
    for(ivar in 1:length(clim_var)){
      vars <- grep(clim_var[ivar],names(clim_data_obs)) # indices for ppt columns
      tmp1 <- clim_alt*colMeans(clim_data_obs) # vector of ppt additions
      # matrix of additions to match DIMS of observed climate
      tmp1 <- matrix(tmp1,NROW(clim_data_obs),NCOL(clim_data_obs),byrow=T) 
      if(colMeans(tmp1)[1]!=tmp1[1,1]) { stop("climate change matrix mis-aligned") }
      clim_data[,vars] <- clim_data_obs[,vars]+tmp1[,vars] # applies ppt additions
    }
    # Make sure climate perturbations have been applied
    if(mean(clim_data$pptLag)==mean(clim_data_obs$pptLag)) {
      stop("climate perturbation not applied")
    }
  }
  
  # Retriev spp-specific scalers
  clim_scalers <- subset(qbm_clim_scalers, yearout==1 & species==doSpp)
  scalers <- clim_scalers[,c("means","sds")]
  
  scaled_climate <- t(apply(X = clim_data[,clim_covariates], MARGIN = 1, 
                          FUN = function(x){(x-scalers$means)/scalers$sds}))
  test <- (clim_data[1,clim_covariates] - scalers$means) / scalers$sds
  test_diff <- sum(test-scaled_climate[1,])
  if(test_diff!=0){ stop("error in climate scaling") }
  
  return(scaled_climate)
}



####
####  RUN CLIMATE CHANGE SIMULATIONS
####

###
### Simulation 1: No Climate Change ----
###
clim_var <- NA
for(ispp in 1:length(spp_list)){
  do_species <- spp_list[ispp]
  params <- fetch_qbm_params(do_species = do_species, 
                             dir_to_find = path2mcmcs, 
                             mean_params = use_mean_params) # returns list of mean params
                                                            # or list of MCMC of params
  current_clim <- perturb_climate(clim_alt, clim_var, clim_data_obs, do_spp = ispp)
  
  # Simulate to equilibrium
  outD <- data.frame(cover=NA, species=NA, year=NA)
  cover <- numeric(tsims)
  cover[1] <- 0.01
  pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
  for(t in 2:tsims){
    if(use_mean_params==TRUE){
      inttmp <- params[["intercept"]]$value
      slopetmp <- params[["coveff"]]$value
      tmpclim <- params[["climeff"]]$value[1:7]
      tmptau <- params[["tau"]]$value
    }
    
    if(use_mean_params==FALSE){
      climeff <- params[["climeff"]]
      coveff <- params[["coveff"]]
      intercept <- params[["intercept"]]
      tau <- params[["tau"]]
      randchain <- sample(x = unique(climeff$Chain), size = 1) # sample random chain
      randiter <- sample(x = unique(climeff$Iteration), size = 1) # sample random MCMC iteration
      inttmp <- subset(intercept, Chain==randchain & Iteration==randiter)$value
      slopetmp <- subset(coveff, Chain==randchain & Iteration==randiter)$value
      tmpclim <- subset(climeff, Chain==randchain & Iteration==randiter)$value[1:7]
      tmptau <- subset(tau, Chain==randchain & Iteration==randiter)$value
    }
    
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- current_clim[climyear,]
    
    cover[t] <- growFunc(N = cover[t-1], int = inttmp, 
                         slope = slopetmp, clims = tmpclim,
                         climcovs = climcovs, tau = tmptau) 
    
    setTxtProgressBar(pb, t)
  }#end simulation loop
  
  # Save the output
  covd <- as.data.frame(cover[burn.in:tsims])
  covd$species <- do_species
  covd$climsim <- "obs"
  colnames(covd)[1] <- "cover"
  
  # Set output filename
  filenameid <- paste(clim_var, collapse="")
  outfile <- paste0(resultspath,spp_list[ispp],"_qbm_cover_", filenameid, "Change.csv")
  write.csv(covd, outfile)
  cat("\n") # makes a space before next line printing
  cat(paste("DONE WITH", spp_list[ispp], "FOR OBSERVED CLIMATE SIMULATION"))
  cat("\n") # makes a space before next line printing
} # end species loop




###
### Simulation 2: Increase precipitation ----
###
clim_var <- "ppt"
for(ispp in 1:length(spp_list)){
  do_species <- spp_list[ispp]
  params <- fetch_qbm_params(do_species = do_species, 
                             dir_to_find = path2mcmcs, 
                             mean_params = use_mean_params) # returns list of mean params
  # or list of MCMC of params
  current_clim <- perturb_climate(clim_alt, clim_var, clim_data_obs, do_spp = ispp)
  
  # Simulate to equilibrium
  outD <- data.frame(cover=NA, species=NA, year=NA)
  cover <- numeric(tsims)
  cover[1] <- 0.01
  pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
  for(t in 2:tsims){
    if(use_mean_params==TRUE){
      inttmp <- params[["intercept"]]$value
      slopetmp <- params[["coveff"]]$value
      tmpclim <- params[["climeff"]]$value[1:7]
      tmptau <- params[["tau"]]$value
    }
    
    if(use_mean_params==FALSE){
      climeff <- params[["climeff"]]
      coveff <- params[["coveff"]]
      intercept <- params[["intercept"]]
      tau <- params[["tau"]]
      randchain <- sample(x = unique(climeff$Chain), size = 1) # sample random chain
      randiter <- sample(x = unique(climeff$Iteration), size = 1) # sample random MCMC iteration
      inttmp <- subset(intercept, Chain==randchain & Iteration==randiter)$value
      slopetmp <- subset(coveff, Chain==randchain & Iteration==randiter)$value
      tmpclim <- subset(climeff, Chain==randchain & Iteration==randiter)$value[1:7]
      tmptau <- subset(tau, Chain==randchain & Iteration==randiter)$value
    }
    
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- current_clim[climyear,]
    
    cover[t] <- growFunc(N = cover[t-1], int = inttmp, 
                         slope = slopetmp, clims = tmpclim,
                         climcovs = climcovs, tau = tmptau) 
    
    setTxtProgressBar(pb, t)
  }#end simulation loop
  
  # Save the output
  covd <- as.data.frame(cover[burn.in:tsims])
  covd$species <- do_species
  covd$climsim <- "obs"
  colnames(covd)[1] <- "cover"
  
  # Set output filename
  filenameid <- paste(clim_var, collapse="")
  outfile <- paste0(resultspath,spp_list[ispp],"_qbm_cover_", filenameid, "Change.csv")
  write.csv(covd, outfile)
  cat("\n") # makes a space before next line printing
  cat(paste("DONE WITH", spp_list[ispp], "FOR +PPT SIMULATION"))
  cat("\n") # makes a space before next line printing
} # end species loop



###
### Simulation 3: Increase temperature ----
###
clim_var <- "Tmean"
for(ispp in 1:length(spp_list)){
  do_species <- spp_list[ispp]
  params <- fetch_qbm_params(do_species = do_species, 
                             dir_to_find = path2mcmcs, 
                             mean_params = use_mean_params) # returns list of mean params
  # or list of MCMC of params
  current_clim <- perturb_climate(clim_alt, clim_var, clim_data_obs, do_spp = ispp)
  
  # Simulate to equilibrium
  outD <- data.frame(cover=NA, species=NA, year=NA)
  cover <- numeric(tsims)
  cover[1] <- 0.01
  pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
  for(t in 2:tsims){
    if(use_mean_params==TRUE){
      inttmp <- params[["intercept"]]$value
      slopetmp <- params[["coveff"]]$value
      tmpclim <- params[["climeff"]]$value[1:7]
      tmptau <- params[["tau"]]$value
    }
    
    if(use_mean_params==FALSE){
      climeff <- params[["climeff"]]
      coveff <- params[["coveff"]]
      intercept <- params[["intercept"]]
      tau <- params[["tau"]]
      randchain <- sample(x = unique(climeff$Chain), size = 1) # sample random chain
      randiter <- sample(x = unique(climeff$Iteration), size = 1) # sample random MCMC iteration
      inttmp <- subset(intercept, Chain==randchain & Iteration==randiter)$value
      slopetmp <- subset(coveff, Chain==randchain & Iteration==randiter)$value
      tmpclim <- subset(climeff, Chain==randchain & Iteration==randiter)$value[1:7]
      tmptau <- subset(tau, Chain==randchain & Iteration==randiter)$value
    }
    
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- current_clim[climyear,]
    
    cover[t] <- growFunc(N = cover[t-1], int = inttmp, 
                         slope = slopetmp, clims = tmpclim,
                         climcovs = climcovs, tau = tmptau) 
    
    setTxtProgressBar(pb, t)
  }#end simulation loop
  
  # Save the output
  covd <- as.data.frame(cover[burn.in:tsims])
  covd$species <- do_species
  covd$climsim <- "obs"
  colnames(covd)[1] <- "cover"
  
  # Set output filename
  filenameid <- paste(clim_var, collapse="")
  outfile <- paste0(resultspath,spp_list[ispp],"_qbm_cover_", filenameid, "Change.csv")
  write.csv(covd, outfile)
  cat("\n") # makes a space before next line printing
  cat(paste("DONE WITH", spp_list[ispp], "FOR +TEMP SIMULATION"))
  cat("\n") # makes a space before next line printing
} # end species loop



###
### Simulation 4: Increase precipitation and temperature ----
###
clim_var <- c("ppt","Tmean")
for(ispp in 1:length(spp_list)){
  do_species <- spp_list[ispp]
  params <- fetch_qbm_params(do_species = do_species, 
                             dir_to_find = path2mcmcs, 
                             mean_params = use_mean_params) # returns list of mean params
  # or list of MCMC of params
  current_clim <- perturb_climate(clim_alt, clim_var, clim_data_obs, do_spp = ispp)
  
  # Simulate to equilibrium
  outD <- data.frame(cover=NA, species=NA, year=NA)
  cover <- numeric(tsims)
  cover[1] <- 0.01
  pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
  for(t in 2:tsims){
    if(use_mean_params==TRUE){
      inttmp <- params[["intercept"]]$value
      slopetmp <- params[["coveff"]]$value
      tmpclim <- params[["climeff"]]$value[1:7]
      tmptau <- params[["tau"]]$value
    }
    
    if(use_mean_params==FALSE){
      climeff <- params[["climeff"]]
      coveff <- params[["coveff"]]
      intercept <- params[["intercept"]]
      tau <- params[["tau"]]
      randchain <- sample(x = unique(climeff$Chain), size = 1) # sample random chain
      randiter <- sample(x = unique(climeff$Iteration), size = 1) # sample random MCMC iteration
      inttmp <- subset(intercept, Chain==randchain & Iteration==randiter)$value
      slopetmp <- subset(coveff, Chain==randchain & Iteration==randiter)$value
      tmpclim <- subset(climeff, Chain==randchain & Iteration==randiter)$value[1:7]
      tmptau <- subset(tau, Chain==randchain & Iteration==randiter)$value
    }
    
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- current_clim[climyear,]
    
    cover[t] <- growFunc(N = cover[t-1], int = inttmp, 
                         slope = slopetmp, clims = tmpclim,
                         climcovs = climcovs, tau = tmptau) 
    
    setTxtProgressBar(pb, t)
  }#end simulation loop
  
  # Save the output
  covd <- as.data.frame(cover[burn.in:tsims])
  covd$species <- do_species
  covd$climsim <- "obs"
  colnames(covd)[1] <- "cover"
  
  # Set output filename
  filenameid <- paste(clim_var, collapse="")
  outfile <- paste0(resultspath,spp_list[ispp],"_qbm_cover_", filenameid, "Change.csv")
  write.csv(covd, outfile)
  cat("\n") # makes a space before next line printing
  cat(paste("DONE WITH", spp_list[ispp], "FOR +PPT +TEMP SIMULATION"))
  cat("\n") # makes a space before next line printing
} # end species loop
