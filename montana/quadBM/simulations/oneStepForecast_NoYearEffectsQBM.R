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
pGrowAll <- subset(pGrow, Coef=="intG"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2"|Coef=="tau"|Coef=="betaSpp"|Coef=="intercept")
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
growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim, doGrp){
  growNow <- subset(pGrowAll, Spp==sppSim)
  iID <- which(growNow$Coef=="intercept")
  intercept <- growNow$value[iID]
  sID <- which(growNow$Coef=="betaSpp")
  size <- growNow$value[sID]
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  tID <- which(growNow$Coef=="tau")
  tau <- growNow$value[tID]
  gID <- which(growNow$Coef=="intG")
  intG <- growNow$value[gID][doGrp]
  newN <- intercept+intG+size*N+sum(climEffs*climate)
  newN <- antilogit(newN)
  print(newN)
  p <- newN * tau
  q <- (1 - newN) * tau
  C <- rbeta(1, p, q)
  print(C)
  return(C)
}


####
#### Run simulations -----------------------------------------------------
####
outDraw <- data.frame(year=NA, cover=NA, sim=NA, species=NA, quad=NA)
for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  nSim <- NumberSimsPerYear
  yearsN <- length(unique(allD$year))
  years <- as.numeric(unique(allD$year)+1900)
  yearsID <- unique(allD$year)
  Nsave <- matrix(ncol=yearsN, nrow=nSim)
  sppD <- subset(allD, Species==sppSim)
  quadList <- as.data.frame(as.character(unique(sppD$quad)))
  quadList$group <- substring(quadList[,1], 1, 1)
  quadList$groupNum <- as.numeric(as.factor(quadList$group))
  colnames(quadList)[1] <- "quad"
  
  for(qd in 1:nrow(quadList)){
    for(yr in 2:yearsN){
      Nstart <- subset(sppD, year==yearsID[yr-1] & quad==quadList[qd,1])$percCover
      for(sim in 1:nSim){
        climate <- subset(climD, year==years[yr-1])[,c(3,5,4,6)]
        Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=Nstart, climate=climate, simsPerYear=length(NforG), doYear=years[yr], sppSim=sppSim, doGrp=quadList[qd,3])
        Nsave[sim,yr] <- Nout
        print(paste("simulation", sim, "of year", yr, "in quad", quadList[qd,1], "for", sppList[i]))
      }#end simulations loop
    }#end year loop
    dN <- as.data.frame(Nsave)
    colnames(dN) <- years
    nM <- melt(dN)
    nM$sim <- rep(1:nSim, length(years))
    nM$species <- rep(sppSim, nSim*length(years))
    nM$quad <- quadList[qd,1]
    colnames(nM)[1:2] <- c("year", "cover")
    outDraw <- rbind(outDraw, nM)
  }#end group loop
}#end species loop

####
#### Make plots
####
outD <- outDraw[2:nrow(outDraw),]
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
allD$year <- allD$year+1900
combD <- merge(outD, allD, by.x = c("species", "year", "quad"), by.y = c("Species", "year", "quad"))
lagD <- combD
lagD$year2 <- as.numeric(lagD$year)+1
lagD <- lagD[,c(1,3,5,8,9)]
colnames(lagD)[4:5] <- c("lagCover", "year")
combD2 <- merge(combD, lagD, by=c("species", "quad", "year", "sim"))
combD2$coverChangeObs <- with(combD2, percCover*100-lagCover*100)
combD2$coverChangePred <- with(combD2, cover*100-lagCover*100)
combD2$resids <- with(combD2, coverChangeObs - coverChangePred)
resD <- combD2

library(ggplot2)
ggplot(data=resD, aes(x=year, y=resids))+
  geom_boxplot()+
  facet_grid(species~., scales = "free")
