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

NumberSimsPerYear <- 10

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
pGrow <- read.csv("../vitalRateRegressions/betaLikelihood/growth/growthStats.csv")
pGrow$Coef <- c(rep("beta", 13*4),
                rep("betaSpp", 4),
                rep("intG", 6*4),
                rep("intYr", 13*4),
                rep("intercept", 4),
                rep("rain1", 4),
                rep("rain2", 4),
                rep("tau", 4),
                rep("temp1", 4),
                rep("temp2", 4))
pGrow <- pGrow[,c(2,6,7)]
colnames(pGrow) <- c("value", "Spp", "Coef")
pGrowAll <- subset(pGrow, Coef=="intG"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2"|Coef=="tau")
pGrowYrs <- subset(pGrow, Coef=="beta" | Coef=="intYr")
years <- unique(allD$year)[2:14]+1900
pGrowYrs$Year <- c(rep(years, each=4),
                   rep(years, each=4))

####
#### Utility functions
####
#antilogit function
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
logit <- function(x) { log(x / (1-x) )}

####
#### Vital rate functions -----------------------------------------------
####
growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim){
  growNow <- subset(pGrowAll, Spp==sppSim)
  growNowYr <- subset(pGrowYrs, Year==doYear)
  growNowYr <- subset(growNowYr, Spp==sppSim)
  iID <- which(growNowYr$Coef=="intYr")
  intercept <- growNowYr$value[iID]
  sID <- which(growNowYr$Coef=="beta")
  size <- growNowYr$value[sID]
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  tID <- which(growNow$Coef=="tau")
  tau <- growNow$value[tID]
  newN <- intercept+size*N+sum(climEffs*climate)
  newN <- antilogit(newN)
#   print(newN)
  p <- newN * tau
  q <- (1 - newN) * tau
  C <- rbeta(1, p, q)
#   print(C)
  return(newN)
}


####
#### Run simulations -----------------------------------------------------
####
outD <- data.frame(variable=NA, cover=NA, sim=NA, species=NA)
length(sppList)
for(i in 1:1){
  sppSim <- sppList[i]
  nSim <- 50
  yearsN <- 100
  years <- unique(allD$year)+1900
  yearsID <- unique(allD$year)
  Nsave <- matrix(ncol=yearsN, nrow=nSim)
  sppD <- subset(allD, Species==sppSim)
  Nsave[,1] <- mean(subset(sppD, year==yearsID[1])$percCover)
  
  
  for(yr in 2:yearsN){
    for(sim in 1:nSim){
      N <- Nsave[sim,yr-1]
      climYr <- sample(climD$year,1)
      climate <- subset(climD, year==climYr)[,c(3,5,4,6)]
      doYear <- sample(years[2:length(years)], 1)
      Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=N, climate=climate, simsPerYear=length(NforG), doYear=doYear, sppSim=sppSim)
      Nsave[sim,yr] <- Nout
      print(paste("Simulation", sim, "of year", yr, "for", sppSim))
    }#end sim loop
  }#end year loop
  
  dN <- as.data.frame(Nsave)
  colnames(dN) <- seq(1:yearsN)
  nM <- melt(dN)
  nM$sim <- rep(1:nSim, length(yearsN))
  nM$species <- rep(sppSim, nSim*length(yearsN))
  colnames(nM)[2] <-  "cover"
  outD <- rbind(outD, nM)
}

####
#### Output
####
outD <- outD[2:nrow(outD),]
outP <- ddply(outD, .(species, as.numeric(variable)), summarize,
              coverAvg = median(cover),
              up = quantile(cover, 0.75))
colnames(outP)[2] <- "year"
plot(outP$year, outP$coverAvg*100, type="l")
