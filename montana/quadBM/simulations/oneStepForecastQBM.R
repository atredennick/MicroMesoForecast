#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

NumberSimsPerYear <- 100

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

doSpp <- sppList[1]

#load vital rate parameters
pCol <- readRDS("../vitalRateRegressions/colonization/colonizationParamsMCMC.rds")
pGrow <- readRDS("../vitalRateRegressions/growth/growthParamsMCMC.rds")
pSurv <- readRDS("../vitalRateRegressions/survival/survivalParamsMCMC.rds")

####
#### Organize parameter values
####
#colonization
pCol2 <- melt(pCol)
pCol2$Spp <- c(rep(rep(sppList, each=3000), times=6),
               rep(rep(sppList, each=3000), times=5))
pCol2$Coef <- c(rep("gInt", times=6*4*3000),
                rep("int", times=4*3000),
                rep("rain1", times=4*3000),
                rep("rain2", times=4*3000),
                rep("temp1", times=4*3000),
                rep("temp2", times=4*3000))
colnames(pCol2)[1] <- "Iter"
pCol <- pCol2[,c(1,3:5)]; rm(pCol2)

#survival
pSurv2 <- melt(pSurv)
pSurv2$Spp <- c(rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=6),
               rep(rep(sppList, each=3000), times=5))
pSurv2$Coef <- c(rep("beta", times=3000*4),
                 rep("gInt", times=6*4*3000),
                rep("int", times=4*3000),
                rep("rain1", times=4*3000),
                rep("rain2", times=4*3000),
                rep("temp1", times=4*3000),
                rep("temp2", times=4*3000))
colnames(pSurv2)[1] <- "Iter"
pSurv <- pSurv2[,c(1,3:5)]; rm(pSurv2)

#growth
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=4))
pGrow2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(pGrow2)[1] <- "Iter"
pGrow <- pGrow2[,c(1,3:5)]; rm(pGrow2)

pGrowAll <- subset(pGrow, Coef=="gInt"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2")
pGrowYrs <- subset(pGrow, Coef=="beta" | Coef=="intYr")
years <- unique(allD$year)[2:14]+1900
pGrowYrs$Year <- c(rep(rep(years, each=3000), each=4),
                   rep(rep(years, each=3000), each=4))



####
#### Utility functions
####
#antilogit function
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
logit <- function(x) { log(x / (1-x) )}

####
#### Vital rate functions -----------------------------------------------
####
survFunc <- function(pSurv, N, climate, simsPerYear, doYear, sppSim){
  survNow <- subset(pSurv, Spp==sppSim)
  doNow <- sample(x = c(1:3000), 1)
  survNow <- subset(survNow, Iter == doNow)
  iID <- which(survNow$Coef=="int")
  intercept <- survNow$value[iID]
  sID <- which(survNow$Coef=="beta")
  size <- survNow$value[sID]
  cID <- which(survNow$Coef=="rain1"|survNow$Coef=="rain2"|survNow$Coef=="temp1"|survNow$Coef=="temp2")
  climEffs <- survNow$value[cID]
  newN <- intercept+size*N+sum(climEffs*climate)
  newN <- antilogit(newN)
  return(newN)
}

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
  newN <- intercept+size*N+sum(climEffs*climate)
  newN <- antilogit(newN)
  return(newN)
}

colFunc <- function(pCol, N, climate, simsPerYear, doYear, sppSim){
  colNow <- subset(pCol, Spp==sppSim)
  doNow <- sample(x = c(1:3000), 1)
  colNow <- subset(colNow, Iter==doNow)
  iID <- which(colNow$Coef=="int")
  intercept <- colNow$value[iID]
  cID <- which(colNow$Coef=="rain1"|colNow$Coef=="rain2"|colNow$Coef=="temp1"|colNow$Coef=="temp2")
  climEffs <- colNow$value[cID]
  newN <- intercept+sum(climEffs*climate)
  newN <- antilogit(newN)
  return(newN*0.01)
}


####
#### Run simulations -----------------------------------------------------
####
outD <- data.frame(year=NA, cover=NA, sim=NA, species=NA)

for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  nSim <- NumberSimsPerYear
  yearsN <- length(unique(allD$year))
  years <- unique(allD$year)+1900
  yearsID <- unique(allD$year)
  Nsave <- matrix(ncol=yearsN, nrow=nSim)
  sppD <- subset(allD, Species==sppSim)
  
  for(yr in 2:yearsN){
    Nstart <- mean(subset(sppD, year==yearsID[yr-1])$percCover)
    for(sim in 1:nSim){
      climate <- subset(climD, year==years[yr-1])[,c(3,5,4,6)]
      
      ifelse(Nstart[Nstart>0],
             Nout <- survFunc(pSurv=pSurv, N=Nstart, climate=climate, simsPerYear=length(NforG), doYear=years[yr], sppSim=sppSim)*growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=Nstart, climate=climate, simsPerYear=length(NforG), doYear=years[yr], sppSim=sppSim),
             Nout <- colFunc(pCol=pCol, N=Nstart, climate=climate, simsPerYear=length(NforC), doYear=years[yr], sppSim=sppSim))
      Nsave[sim,yr] <- Nout
      print(paste("simulation", sim, "of year", yr, "for", sppList[i]))
    }#end simulations loop
  }#end year loop
  
  dN <- as.data.frame(Nsave)
  colnames(dN) <- years
  nM <- melt(dN)
  nM$sim <- rep(1:nSim, length(years))
  nM$species <- rep(sppSim, nSim*length(years))
  colnames(nM)[1:2] <- c("year", "cover")
  outD <- rbind(outD, nM)
}#end species loop

####
#### Make plots
####
outD <- outD[2:nrow(outD),]
quadD1 <- ddply(allD, .variables = c("year", "Species"), .fun = summarise,
               year = mean(year),
               cover = mean(percCover))
sppD <- ddply(allD, .variables = c("Species"), .fun = summarise,
              avgcover = mean(percCover))
quadD <- merge(quadD1, sppD)
quadD$year <- quadD$year+1900
quadD$yearDiff <- with(quadD, cover*100-avgcover*100)
tmp <- which(outD$year==1932)
tmp2 <- which(quadD$year==1932)
outD[tmp, "cover"] <- quadD[tmp2, "cover"]

#create lag cover variable
lagD <- outD
lagD$year2 <- as.numeric(lagD$year)+1
lagD <- lagD[, c(2,4,5)]
colnames(lagD) <- c("lagCover", "species", "year")
outD2 <- merge(outD, lagD, by=c("species", "year"))

lagQ <- quadD
lagQ$year2 <- as.numeric(lagQ$year)+1
lagQ <- lagQ[, c(1,3,6)]
colnames(lagQ) <- c( "Species", "lagCover", "year")
quadD2 <- merge(quadD, lagQ, by=c("Species", "year"))

#calculate expected and observed cover change
quadD2$coverChange <- with(quadD2, cover*100-lagCover*100)
outD2$coverChange <- with(outD2, cover*100-lagCover*100)

resD <- merge(outD2, quadD2, by.x = c("species", "year"), by.y = c("Species", "year"))
resD$covChangeResiduals <- with(resD, coverChange.x - lagCover.y)
resD <- subset(resD, year!=1932)
saveRDS(resD, file = "quadBM_oneStep_Residuals.rds")

