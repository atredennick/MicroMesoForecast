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
NumberSimsPerYear <- 100

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

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
growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim, doGrp){
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
outDraw <- data.frame(quad=NA, cover=NA, sim=NA, species=NA, year=NA)
for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  nSim <- NumberSimsPerYear
  yearsN <- length(unique(allD$year))
  years <- as.numeric(unique(allD$year)+1900)
  yearsID <- unique(allD$year)
  sppD <- subset(allD, Species==sppSim)
  
  for(yr in 2:yearsN){
    yrD <- subset(sppD, year==yearsID[yr-1])
    quadList <- as.data.frame(as.character(unique(yrD$quad)))
    quadList$group <- substring(quadList[,1], 1, 1)
    quadList$groupNum <- as.numeric(as.factor(quadList$group))
    colnames(quadList)[1] <- "quad"
    climate <- subset(climD, year==years[yr-1])[,c(3,5,4,6)]
    Nsave <- matrix(ncol=nrow(quadList), nrow=nSim)
    for(qd in 1:nrow(quadList)){
      Nstart <- subset(yrD, quad==as.character(quadList[qd,1]))$percCover
      for(sim in 1:nSim){
        Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=Nstart, climate=climate, doYear=years[yr], sppSim=sppSim, doGrp=quadList[qd,3])
        Nsave[sim,qd] <- Nout
        print(paste("simulation", sim, "of year", yr, "in quad", quadList[qd,1], "for", sppList[i]))
      }#end simulations loop
    }#end group loop
    dN <- as.data.frame(Nsave)
    colnames(dN) <- as.character(quadList[,1])
    nM <- melt(dN)
    nM$sim <- rep(1:nSim, nrow(quadList))
    nM$species <- rep(sppSim, nSim*nrow(quadList))
    nM$year <- years[yr]
    colnames(nM)[1:2] <- c("quad", "cover")
    outDraw <- rbind(outDraw, nM)
  }#end year loop
}#end species loop

####
#### Make plots
####
outD <- outDraw[2:nrow(outDraw),]
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
allD$year <- allD$year+1900
combD <- merge(outD, allD, by.x = c("species", "year", "quad"), by.y = c("Species", "year", "quad"), all.y=TRUE)
id <- which(is.na(combD$sim)==TRUE)
combD$sim[id] <- 1
lagD <- subset(combD, sim==1)
lagD$year2 <- as.numeric(lagD$year)+1
lagD <- lagD[,c(1,3,8,9)]
colnames(lagD)[3:4] <- c("lagCover", "year")
combD2 <- merge(combD, lagD, by=c("species", "quad", "year"))
combD2$coverChangeObs <- with(combD2, percCover*100-lagCover*100)
combD2$coverChangePred <- with(combD2, cover*100-lagCover*100)
combD2$resids <- with(combD2, coverChangeObs - coverChangePred)
resD <- combD2
saveRDS(resD, file = "oneStepResids_YearEffects.rds")
# noZ <- which(resD$lagCover!=0)
# resD <- resD[noZ,]
# 
# library(ggplot2)
# ggplot(data=resD, aes(x=year, y=resids, group=year))+
#   geom_boxplot()+
#   facet_grid(species~., scales = "free")
