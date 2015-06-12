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
rm(list=ls(all=TRUE))

## Set do_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_year <- as.numeric(myargument)

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
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

##  Climate
climD <- read.csv("../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)


####
####  Read in statistical model parameters ------------
####
in_file <- paste("../vitalRateRegressions/truncNormModel/validation/results/growthParamsMCMC_", do_year, ".rds", sep="")
pGrow <- readRDS(in_file)
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=12),
                rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=12),
                rep(rep(sppList, each=3000), times=7))
pGrow2$Coef <- c(rep("beta", times=3000*4*12),
                 rep("betaSpp", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*12),
                 rep("intercept", times=3000*4),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("rainLag", times=4*3000),
                 rep("tau", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(pGrow2)[1] <- "Iter"
pGrow <- pGrow2[,c(1,3:5)]; rm(pGrow2)
pGrowAll <- subset(pGrow, Coef=="gInt"|Coef=="rain1"|
                          Coef=="rain2"|Coef=="rainLag"|
                          Coef=="temp1"|Coef=="temp2"|
                          Coef=="tau"|Coef=="betaSpp"|Coef=="intercept")


####
#### Population growth function ----------------------------
####
growFunc <- function(pGrowAll, N, climate, 
                     simsPerYear, sppSim, doGrp){
  growNow <- subset(pGrowAll, Spp==sppSim)
  doNow <- sample(x = c(1:3000), 1)
  growNow <- subset(growNow, Iter==doNow)
  iID <- which(growNow$Coef=="intercept")
  intercept <- growNow$value[iID]
  sID <- which(growNow$Coef=="betaSpp")
  size <- growNow$value[sID]
  cID <- which(growNow$Coef=="rainLag"|growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  gID <- which(growNow$Coef=="intG")
  intG <- growNow$value[gID][doGrp]
  tID <- which(growNow$Coef=="tau")
  tau <- growNow$value[tID]
  mu <- intercept+size*log(N)+sum(climEffs*climate)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = sqrt(1/tau), min = 0, max = 1)
  return(newN)
}


####
#### Run simulations -----------------------------------------------------
####
years <- as.numeric(unique(allD$year)+1900)
yearsID <- unique(allD$year)
yrD <- subset(allD, year==(do_year-1))
climate <- subset(climD, year==(1900+do_year-1))[,c(2,3,5,4,6)]
nSim <- NumberSimsPerYear

outDraw <- data.frame(quad=NA, cover=NA, sim=NA, species=NA, year=NA)
for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  sppD <- subset(yrD, Species==sppSim)
  quadList <- as.data.frame(as.character(unique(sppD$quad)))
  quadList$group <- substring(quadList[,1], 1, 1)
  quadList$groupNum <- as.numeric(as.factor(quadList$group))
  colnames(quadList)[1] <- "quad"
  
  Nsave <- matrix(ncol=nrow(quadList), nrow=nSim)
  for(qd in 1:nrow(quadList)){
    Nstart <- subset(sppD, quad==as.character(quadList[qd,1]))$percCover
    if(Nstart>0){
      for(sim in 1:nSim){
        Nout <- growFunc(pGrow=pGrowAll, N=Nstart, climate=climate, sppSim=sppSim, doGrp=quadList[qd,3])
        Nsave[sim,qd] <- Nout
        #       print(paste("simulation", sim, "of year", yr, "in quad", quadList[qd,1], "for", sppList[i]))
      }#end simulations loop
    }#End empty quad if statement
  }#end group loop
  dN <- as.data.frame(Nsave)
  colnames(dN) <- as.character(quadList[,1])
  nM <- melt(dN)
  nM$sim <- rep(1:nSim, nrow(quadList))
  nM$species <- rep(sppSim, nSim*nrow(quadList))
  nM$year <- do_year
  colnames(nM)[1:2] <- c("quad", "cover")
  outDraw <- rbind(outDraw, nM)
}#end species loop

sims_to_remove <- which(is.na(outDraw$cover)==TRUE)
loyo_simulations <- outDraw[-sims_to_remove,]
colnames(loyo_simulations) <- c("quad", "predicted_cover",
                                "sim", "species", "predicted_year")

####
####  Save output ---------------------------------
####
out_file <- paste("loyo_validation_sims_", do_year, ".rds", sep="")
saveRDS(loyo_simulations, out_file)
