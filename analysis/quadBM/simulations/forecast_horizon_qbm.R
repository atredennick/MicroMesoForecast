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
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../weather/Climate.csv")

backD <- data.frame(climYear=NA,
                    quad = NA,
                    year= NA,
                    totCover= NA,
                    Species= NA,
                    propCover= NA,
                    lag.cover= NA,
                    pptLag= NA,
                    ppt1= NA,
                    TmeanSpr1= NA,
                    ppt2= NA,
                    TmeanSpr2= NA,
                    TmeanSum1= NA,
                    TmeanSum2= NA,
                    yearID= NA,
                    group = NA,
                    percCover = NA,
                    percLagCover = NA)

#loop through species and remake data frame
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  backD <- rbind(backD, growD)
}#end species loop
obs_data <- backD[2:nrow(backD),]
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
  
  outall <- data.frame(quad=NA, finalcover=NA, sim=NA, 
                       species=NA, yearstart=NA, obs_finalcover=NA)
  for(do_year in all_years){
    ####
    #### Run simulations -----------------------------------------------------
    ####
    yrD <- subset(obs_data, year==do_year & Species==do_species)
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
      Nstart <- subset(yrD, quad==as.character(quadList[qd,1]))$percCover
      if(Nstart>0){
        for(sim in 1:nSim){
          Nnow <- Nstart
          for(yearsim in years2sim){
            print(Nnow)
            weather_year <- weather_scaled[length(all_years)-(max(all_years)-yearsim),]
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
    obs_by_quad <- subset(obs_data, year==45 & Species==do_species)[,c("quad","propCover")]
    colnames(obs_by_quad) <- c("quad", "obs_finalcover")
    tmpoutall <- merge(newNs, obs_by_quad)
    outall <- rbind(outall, tmpoutall)
  }#end year loop
  
  outname <- paste(do_species,"_final_year_cover.RDS", sep="")
  saveRDS(outall, paste("./results/forecast_horizon/", outname, sep=""))
}#end species loop

