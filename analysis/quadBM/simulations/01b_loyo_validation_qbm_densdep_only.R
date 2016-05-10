##  This script runs validation models for the quad-based population model.

##  Statistical results are from leave-one-year-out fits, so here we attempt
##    to use the model results to predict the year left out of fitting. We
##    fit 13 different models, so this script loops through the 13 unfit years
##    to predict to cover value for that year. We also loop through 100 different
##    simulations for each quadrat in each year to capture parameter uncertainty 
##    by drawing parameter values from the MCMC chain for each simulation. For
##    the simulation we do not use the random year effects -- those are just for
##    making sure our climate coefficients are unbiased.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

##  Date began:     3.23.2015
##  Date completed: 3.23.2015
##  Date tested:    3.24.2015 -- Non-exhaustive tests, just made sure 'if' statement
##                                is working and that the growth function is
##                                operating as expected.

##  The script takes 'do_year' as a command line prompt. So, e.g.,
##    run as: "R CMD BATCH -33 loyo_validation_qbm.R" for year 33.

##  Depends: script needs the following R packages
##            -- install.packages(c("reshape2", "plyr", "EnvStats"))


##  Clear the workspace
rm(list=ls(all.names=TRUE))

## Set do_year for validation from command line prompt
# args <- commandArgs(trailingOnly = F)
# myargument <- args[length(args)]
# myargument <- sub("-","",myargument)
# do_year <- as.numeric(myargument)

##  Set number of simulations per year
NumberSimsPerYear <- 100

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
obs_data <- readRDS("../../processed_data/cover_with_weather.RDS")
all_years <- unique(obs_data$year)
sppList <- unique(obs_data$species)

####
#### Population growth function ----------------------------
####
##  Source cover change, population growth function
source("qbm_sim_fxn.R")


##  Start looping over species and year within species
for(do_species in sppList){
  outall <- list()
  for(do_year in all_years){
    ####
    ####  Read in statistical model parameters ------------
    ####
    yearnow <- do_year
    fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/validation/fits_noclimate/popgrowth_stanmcmc_noclimate_", 
                             do_species, "_leaveout", yearnow, ".RDS", sep=""))
    ##  Break up MCMC into regression components
    fitthin <- fitlong
    
    # Mean cover (size) effects
    coveff <- fitthin[grep(glob2rx("b1_mu"), fitthin$Parameter),]
    
    # Mean intercepts
    intercept <- fitthin[grep("a_mu", fitthin$Parameter),]
    
    # Group offsets
    goffs <- fitthin[grep("gint", fitthin$Parameter),]
    goffs$groupid <- substr(goffs$Parameter, 6, length(goffs$Parameter))
    goffs$groupid <- unlist(strsplit(goffs$groupid, split=']'))
    
    # Lognormal sigma (called tau here)
    tau <- fitthin[grep("tau", fitthin$Parameter),]
    
    
    ####
    #### Run simulations -----------------------------------------------------
    ####
    yrD <- subset(obs_data, year==yearnow & species==do_species)
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
          randchain <- sample(x = coveff$Chain, size = 1)
          randiter <- sample(x = coveff$Iteration, size = 1)
          inttmp <- subset(intercept, Chain==randchain & 
                             Iteration==randiter)
          grptmp <- subset(goffs, Chain==randchain & 
                             Iteration==randiter &
                             groupid==quadList[qd,"groupNum"])
          slopetmp <- subset(coveff, Chain==randchain & 
                               Iteration==randiter)
          tmptau <- subset(tau, Chain==randchain & 
                             Iteration==randiter)
          Nout <- growFunc_noClim(N = Nstart, int = inttmp$value+grptmp$value, 
                           slope = slopetmp$value, tau = tmptau$value) 
          Nsave[sim,qd] <- Nout
          Nstarts[sim,qd] <- Nstart
          print(paste("simulation", sim, "of year", yearnow, "in quad", quadList[qd,1], "for", do_species))
        }#end simulations loop
      }#End empty quad if statement
    }#end group loop
    newNs <- as.data.frame(Nsave)
    colnames(newNs) <- as.character(quadList[,1])
    newNs <- melt(newNs)
    newNs$sim <- rep(1:nSim, nrow(quadList))
    newNs$species <- rep(do_species, nSim*nrow(quadList))
    newNs$year <- yearnow
    colnames(newNs)[1:2] <- c("quad", "pred_cover.t1")
    
    startNs <- as.data.frame(Nstarts)
    colnames(startNs) <- as.character(quadList[,1])
    startNs <- melt(startNs)
    startNs$sim <- rep(1:nSim, nrow(quadList))
    startNs$species <- rep(do_species, nSim*nrow(quadList))
    startNs$year <- yearnow
    colnames(startNs)[1:2] <- c("quad", "obs_cover.t0")
    
    outsaves <- merge(newNs, startNs)
    
    obsNs <- yrD[,c("year", "species", "quad", "propCover.t1")]
    names(obsNs)[which(names(obsNs)=="propCover.t1")] <- "obs_cover.t1"
    outsaves2 <- merge(outsaves, obsNs)
    
    outall <- rbind(outall, outsaves2)
  }#end year loop
  outname <- paste(do_species,"_sim_cover_1step_ahead_year.RDS", sep="")
  saveRDS(outall, paste("./results/one_step_validation_densdep_only/", outname, sep=""))
}#end species loop



