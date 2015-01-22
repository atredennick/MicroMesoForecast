#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(EnvStats)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

#perturb climate data
# climD <- read.csv("../../weather/Climate.csv")
# climScale <- scale(climD[3:6], center = TRUE, scale = TRUE)
# climAvg <- apply(X = climD, MARGIN = 2, FUN = mean)
# climSD <- apply(X = climD, MARGIN = 2, FUN = sd)
# climD[c(3,5)] <- climD[c(3,5)]+(climD[c(3,5)]*0.01)
# climD[c(4,6)] <- climD[c(4,6)]+(climD[c(4,6)]*0.01)
# climD[3] <- (climD[3] - climAvg[3])/climSD[3]
# climD[4] <- (climD[4] - climAvg[4])/climSD[4]
# climD[5] <- (climD[5] - climAvg[5])/climSD[5]
# climD[6] <- (climD[6] - climAvg[6])/climSD[6]

#load vital rate parameters
pGrow <- readRDS("../vitalRateRegressions/truncNormModel/growthParamsMCMC.rds")
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=6))
pGrow2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("betaSpp", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("intercept", times=3000*4),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("tau", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(pGrow2)[1] <- "Iter"
pGrow <- pGrow2[,c(1,3:5)]; rm(pGrow2)

pGrowAll <- subset(pGrow, Coef=="gInt"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2"|Coef=="tau")
pGrowYrs <- subset(pGrow, Coef=="beta" | Coef=="intYr")
years <- unique(allD$year)[2:14]+1900
pGrowYrs$Year <- c(rep(rep(years, each=3000), each=4),
                   rep(rep(years, each=3000), each=4))

####
#### Vital rate functions -----------------------------------------------
####
growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim){
  growNow <- subset(pGrowAll, Spp==sppSim)
  doNow <- sample(x = c(1:3000), 1)
  growNow <- subset(growNow, Iter==doNow)
  growNowYr <- subset(pGrowYrs, Year==doYear)
  growNowYr <- subset(growNowYr, Iter==doNow)
  growNowYr <- subset(growNowYr, Spp==sppSim)
  iID <- which(growNowYr$Coef=="intYr")
  intercept <- growNowYr$value[iID]
  sID <- which(growNowYr$Coef=="beta")
  size <- growNowYr$value[sID]
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  tID <- which(growNow$Coef=="tau")
  tau <- growNow$value[tID]
  mu <- intercept+size*log(N)+sum(climEffs*climate)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = sqrt(1/tau), min = 0, max = 1)
  return(newN)
}

# growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim){
#   growNow <- subset(pGrowAll, Spp==sppSim)
#   growNowYr <- subset(pGrowYrs, Year==doYear)
#   growNowYr <- subset(growNowYr, Spp==sppSim)
#   iID <- which(growNowYr$Coef=="intYr")
#   intercept <- growNowYr$value[iID]
#   sID <- which(growNowYr$Coef=="beta")
#   size <- growNowYr$value[sID]
#   cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
#   climEffs <- growNow$value[cID]
#   tID <- which(growNow$Coef=="tau")
#   tau <- growNow$value[tID]
#   mu <- intercept+size*log(N)+sum(climEffs*climate)
#   newN <- rlnormTrunc(1, meanlog = mu, sdlog = (1/tau), min = 0, max = 1)
# #   newN <- exp(newN)
#   return(newN)
# }


####
#### Run simulations -----------------------------------------------------
####
outD <- data.frame(cover=NA, species=NA, year=NA)

for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  nSim <- 1
  yearsN <- 10000
  years <- unique(allD$year)+1900
  yearsID <- unique(allD$year)
  Nsave <- numeric(yearsN)
  sppD <- subset(allD, Species==sppSim)
  Nsave[1] <- mean(subset(sppD, year==yearsID[1])$percCover)
  
  
  for(yr in 2:yearsN){
    #for(sim in 1:nSim){
      N <- Nsave[yr-1]
      climYr <- sample(climD$year,1)
      climate <- subset(climD, year==climYr)[,c(3,5,4,6)]
      doYear <- sample(years[2:length(years)], 1)
      Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=N, climate=climate, simsPerYear=length(NforG), doYear=doYear, sppSim=sppSim)
      Nsave[yr] <- Nout
      print(paste("Year", yr, "for", sppSim))
    #}#end sim loop
  }#end year loop
  
  dN <- as.data.frame(Nsave)
  colnames(dN) <- "cover"
  dN$species <- rep(sppSim, length(yearsN))
  dN$year <- seq(1:yearsN)
  outD <- rbind(outD, dN)
}

####
#### Output
####
outD <- outD[2:nrow(outD),]
saveRDS(outD, "baselineSimulationQBM.rds")




