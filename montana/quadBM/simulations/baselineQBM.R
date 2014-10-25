#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(reshape2)
library(plyr)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")

doSpp <- sppList[1]

#load vital rate parameters
files <- list.files("../vitalRateRegressions/finalModels/vitalRateParameters/")
paramFiles <- files[grep(doSpp, files)]
pCol <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[1], sep=""))
pGrow <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[2], sep=""))
pSurv <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[3], sep=""))

#antilogit function
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }

####
#### Vital rate functions -----------------------------------------------
####
survFunc <- function(pSurv, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pSurv$fixed[1,1], pSurv$fixed[1,2])
  size <- rnorm(simsPerYear, pSurv$fixed[2,1], pSurv$fixed[2,2])
  yearID <- which(pSurv$yearInt$ID==doYear)
  yearOff <- rnorm(simsPerYear, pSurv$yearInt[yearID,2], pSurv$yearInt[yearID,3])
  yearSize <- rnorm(simsPerYear, pSurv$yearSize[yearID,2], pSurv$yearSize[yearID,3])
  
  newN <- intercept+yearOff+(size+yearSize)*N
  newN <- antilogit(newN)
  return(newN)
}

growFunc <- function(pGrow, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pGrow$fixed[1,1], pGrow$fixed[1,2])
  size <- rnorm(simsPerYear, pGrow$fixed[2,1], pGrow$fixed[2,2])
  yearID <- which(pGrow$yearInt$ID==doYear)
  yearOff <- rnorm(simsPerYear, pGrow$yearInt[yearID,2], pGrow$yearInt[yearID,3])
  yearSize <- rnorm(simsPerYear, pGrow$yearSize[yearID,2], pGrow$yearSize[yearID,3])
  
  newN <- intercept+yearOff+(size+yearSize)*N
  newN <- antilogit(newN)
  return(newN)
}

colFunc <- function(pCol, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pCol$fixed[1,1], pCol$fixed[1,2])
  climCov <- matrix(nrow=simsPerYear, ncol=(nrow(pCol$fixed)-1))
#   for(jj in 1:ncol(climCov)){
#     climCov[,jj] <- rnorm(simsPerYear, pCol$fixed[jj+1,1], pCol$fixed[jj+1,2])
#   } 
  climCov[,1] <- rnorm(simsPerYear, pCol$fixed[2,1], pCol$fixed[2,2])
  colnames(climCov) <- rownames(pCol$fixed)[2:nrow(pCol$fixed)]
  climVars <- numeric(ncol(climCov))
  for(ll in 1:length(climVars)){
    climVars[ll] <- subset(climate, variable==colnames(climCov)[ll])$value
  }
  yearOff <- rnorm(simsPerYear, mean(pCol$yearInt[,2]), mean(pCol$yearInt[,3]))
  
  newN <- intercept+yearOff+sum(climVars*climCov)
  newN <- antilogit(newN)
  return(newN*0.01)
}


####
#### Run simulations -----------------------------------------------------
####
nSim <- 1000
yearsN <- length(unique(allD$year))
years <- unique(allD$year)+1900
yearsID <- unique(allD$year)
Nsave <- matrix(ncol=yearsN, nrow=nSim)
Nsave[,1] <- mean(subset(allD, year==yearsID[1])$propCover)

for(yr in 2:yearsN){
  N <- Nsave[,yr-1]
  climate <- subset(climD, year==years[yr])
  climate <- melt(climate)
  
  NforG <- N[N>0]
  tmpN <- survFunc(pSurv=pSurv, N=NforG, climate=climate, simsPerYear=length(NforG), doYear=yearsID[yr])*growFunc(pGrow=pGrow, N=NforG, climate=climate, simsPerYear=length(NforG), doYear=yearsID[yr])
  
  NforC <- N[N==0]
  colN <- colFunc(pCol=pCol, N=NforC, climate=climate, simsPerYear=length(NforC), doYear=yearsID[yr])
  
  Nout <- c(tmpN, colN)
  Nsave[,yr] <- Nout
}

Nplot <- apply(Nsave, MARGIN = 2, FUN = mean)
plot(years, Nplot*100)
lines(years, Nplot*100)

dN <- as.data.frame(Nsave)
colnames(dN) <- years
nM <- melt(dN)
nM$sim <- rep(1:nSim, length(years))

ggplot(nM, aes(x=variable, y=value*100, group=sim))+
  geom_line(alpha=0.02, color="purple")+
  theme_bw()+
  scale_y_continuous(limits=c(0,25))




