######################################
#### Quad-Based Model simulations ####
######################################

################################################
#### Andrew Tredennick (atredenn@gmail.com) ####
#### 1-7-2014                               ####
################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## This is for one-step-ahead forecasting, by quadrat.                                         ##
## So, initilize with specific quadrat cover, use model with group effect and project forward. ##
## Do this N times for each yearly transition with different MCMC parameters.                  ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#Set number of sims per yearly transition
NumberSimsPerYear <- 100

####
#### Load libraries ----------------------------------
####
library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)

####
#### Bring in data -----------------------------------
####
#cover data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
sppList <- as.character(unique(allD$Species))

#climate data
climD <- read.csv("../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

#load vital rate parameters
pGrow <- readRDS("../vitalRateRegressions/popGrowthParamsMCMC.rds")


####
#### Organize parameter values ---------------------------------
####
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=4),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=4)
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=4))
pGrow2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("betaSpp", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intercept", times=4*3000)
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
#### Vital rate functions -----------------------------------------------
####
growFunc <- function(pGrowAll, pGrowYrs, N, climate, simsPerYear, doYear, sppSim, doGrp){
  growNow <- subset(pGrowAll, Spp==sppSim)
  doNow <- sample(x = c(1:3000), 1)
  growNow <- subset(growNow, Iter==doNow)
  iID <- which(growNow$Coef=="intercept")
  intercept <- growNow$value[iID]
  gID <- which(growNow$Coef=="gInt")
  intG <- growNow$value[gID]
  sID <- which(growNow$Coef=="betaSpp")
  size <- growNow$value[sID]
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  newN <- intercept+intG[doGrp]+size*N+sum(climEffs*climate)
  newN <- exp(newN)
  return(newN)
}


####
#### Run simulations -----------------------------------------------------
####
outDraw <- data.frame(year=NA, cover=NA, sim=NA, species=NA, quad=NA)
length(sppList)
for(i in 1:1){
  sppSim <- sppList[i]
  nSim <- NumberSimsPerYear
  yearsN <- length(unique(allD$year))
  years <- unique(allD$year)+1900
  yearsID <- unique(allD$year)
  Nsave <- matrix(ncol=yearsN, nrow=nSim)
  sppD <- subset(allD, Species==sppSim)
  quadList <- as.data.frame(as.character(unique(sppD$quad)))
  quadList$group <- substring(quadList[,1], 1, 1)
  quadList$groupNum <- as.numeric(as.factor(quadList$group))
  colnames(quadList)[1] <- "quad"
  nrow(quadList)
  for(qd in 1:4){
    for(yr in 2:yearsN){
      Nstart <- round(subset(sppD, year==yearsID[yr-1] & quad==quadList[qd,1])$percCover*100)
      for(sim in 1:nSim){
        climate <- subset(climD, year==years[yr-1])[,c(3,5,4,6)]
        ifelse(length(Nstart)==0,
               Nout <- NA,
               Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=Nstart, climate=climate, simsPerYear=length(NforG), doYear=years[yr], sppSim=sppSim, doGrp=quadList[qd,3]))
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
#### Make output data frame ------------------------------------------------
####
# noNA <- which(is.na(outD$cover)!=TRUE)
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
combD2$coverChangePred <- with(combD2, cover-lagCover*100)
combD2$resids <- with(combD2, coverChangeObs - coverChangePred)
resD <- combD2
saveRDS(resD, file = "oneStepResultsClimateOnly.rds")
# library(ggplot2)
# ggplot(data=resD, aes(x=year, y=resids))+
#   geom_boxplot()+
#   facet_grid(species~., scales = "free")
