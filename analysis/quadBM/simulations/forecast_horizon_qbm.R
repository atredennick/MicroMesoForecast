##  R Script for QPM forecasts of final observation year from all possible
##  lag years. QPM is initialized with cover at year t-x in quadrat q, and 
##  then projected using observed climate to year t. We can then assess the 
##  model's forecast horizon and compare it to the IPM.
##
##  Author: Andrew Tredennick
##  Last update: 4-19-2016


# Clear the workspace
rm(list=ls())

##  Set number of simulations per year
use_mean_parameters <- TRUE
NumberSimsPerYear <- 50
scalers <- readRDS("../../qbm_all_clim_scalers.RDS")

####
####  Load libraries -----------------------------------
####
library('reshape2')
library('plyr')
library('EnvStats')


####
####  Read in data -------------------------------------
####
##  Observations
climD <- read.csv("../../weather/Climate.csv")
obs_data <- readRDS("../../processed_data/cover_with_weather.RDS")
sppList <- sort(unique(obs_data$species))
all_years <- unique(obs_data$year)


####
#### Population growth function ----------------------------
####
##  Source cover change, population growth function
source("qbm_sim_fxn.R")


##  Start looping over species and year within species
for(do_species in sppList){
  ####
  ####  Read in statistical model parameters ------------
  ####
  ##  Scale climate time series
  clim_scalers <- subset(scalers, species==do_species & yearout==1)
  clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
  weather <- climD[,clim_vars]
  weather$ppt1TmeanSpr1 <- weather$ppt1*weather$TmeanSpr1
  weather$ppt2TmeanSpr2 <- weather$ppt2*weather$TmeanSpr2
  weather_scaled <- weather
  for(iclim in 1:nrow(weather)){
    weather_scaled[iclim,] <- (weather[iclim,] - clim_scalers[,"means"])/clim_scalers[,"sds"]
  }
  
  
  fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", 
                           do_species, ".RDS", sep=""))
  ##  Break up MCMC into regression components
  fitthin <- fitlong
  # Climate effects
  climeff <- fitthin[grep("b2", fitthin$Parameter),]
  
  # Mean cover (size) effects
  coveff <- fitthin[grep(glob2rx("b1_mu"), fitthin$Parameter),]
  
  # Mean intercepts
  intercept <- fitthin[grep("a_mu", fitthin$Parameter),]
  
  # Group offsets
  goffs <- fitthin[grep("gint", fitthin$Parameter),]
  goffs$groupid <- substr(goffs$Parameter, 6, length(goffs$Parameter))
  goffs$groupid <- unlist(strsplit(goffs$groupid, split='[.]'))
  
  # Lognormal sigma (called tau here)
  tau <- fitthin[grep("tau", fitthin$Parameter),]
  
  # Set mean parameters if necessary
  if(use_mean_parameters == TRUE) {
    tmpclim <- ddply(climeff, .(Parameter), summarise,
                     avg_value = mean(value))
    names(tmpclim) <- c("Parameter", "value")
    
    inttmp <- ddply(intercept, .(Parameter), summarise,
                    avg_value = mean(value))
    names(inttmp) <- c("Parameter", "value")
    
    grptmp <-  ddply(goffs, .(Parameter), summarise,
                 avg_value = mean(value))
    names(grptmp) <- c("Parameter", "value")
    
    slopetmp <-  ddply(coveff, .(Parameter), summarise,
                     avg_value = mean(value))
    names(slopetmp) <- c("Parameter", "value")
    
    tmptau <-  ddply(tau, .(Parameter), summarise,
                       avg_value = mean(value))
    names(tmptau) <- c("Parameter", "value")
  }
  
  outall <- list()
  for(do_year in all_years){
    ####
    #### Run simulations -----------------------------------------------------
    ####
    yrD <- subset(obs_data, year==do_year & species==do_species)
    years2sim <- do_year:max(all_years)
    nSim <- NumberSimsPerYear
    
    quadList <- as.data.frame(as.character(unique(yrD$quad)))
    quadList$group <- substring(quadList[,1], 1, 1)
    quadList$groupNum <- as.numeric(as.factor(quadList$group))
    colnames(quadList)[1] <- "quad"
    Nsave <- matrix(ncol=nrow(quadList), nrow=nSim)
    Nstarts <- matrix(ncol=nrow(quadList), nrow=nSim)
    
    # Loop over quads to make one-step forecast nSim times
    for(qd in 1:nrow(quadList)){
      Nstart <- subset(yrD, quad==as.character(quadList[qd,1]))$propCover.t0
      if(Nstart>0){
        for(sim in 1:nSim){
          Nnow <- Nstart
          for(yearsim in years2sim){
            # print(Nnow)
            weather_year <- weather_scaled[length(all_years)-(max(all_years)-yearsim),]
            
            if(use_mean_parameters == FALSE) {
              randchain <- sample(x = climeff$Chain, size = 1)
              randiter <- sample(x = climeff$Iteration, size = 1)
              inttmp <- subset(intercept, Chain==randchain & 
                                 Iteration==randiter)
              grptmp <- subset(goffs, Chain==randchain & 
                                 Iteration==randiter &
                                 groupid==quadList[qd,"groupNum"])
              slopetmp <- subset(coveff, Chain==randchain & 
                                   Iteration==randiter)
              tmpclim <- subset(climeff, Chain==randchain & 
                                  Iteration==randiter)
              tmptau <- subset(tau, Chain==randchain & 
                                 Iteration==randiter)
            }
            
            Nnow <- growFunc(N = Nnow, int = inttmp$value+grptmp$value, 
                             slope = slopetmp$value, clims = tmpclim$value,
                             climcovs = weather_year, tau = tmptau$value) 
          }
          Nsave[sim,qd] <- Nnow
        }#end simulations loop
      }#End empty quad if statement
      print(paste("year", do_year, "in quad", quadList[qd,1], "for", do_species))
    }#end group loop
    newNs <- as.data.frame(Nsave)
    colnames(newNs) <- as.character(quadList[,1])
    newNs <- melt(newNs)
    newNs$sim <- rep(1:nSim, nrow(quadList))
    newNs$species <- rep(do_species, nSim*nrow(quadList))
    newNs$yearstart <- do_year
    colnames(newNs)[1:2] <- c("quad", "finalcover")
    obs_by_quad <- subset(obs_data, year==1944 & species==do_species)[,c("quad","propCover.t1")]
    colnames(obs_by_quad) <- c("quad", "obs_finalcover")
    tmpoutall <- merge(newNs, obs_by_quad)
    outall <- rbind(outall, tmpoutall)
  }#end year loop
  
  outname <- paste(do_species,"_final_year_cover.RDS", sep="")
  saveRDS(outall, paste("./results/forecast_horizon/", outname, sep=""))
}#end species loop

