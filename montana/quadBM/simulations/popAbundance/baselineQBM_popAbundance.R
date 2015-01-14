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
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
sppList <- as.character(unique(allD$Species))

#climate data
climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

#load vital rate parameters
pGrow <- readRDS("../../vitalRateRegressions/popAbundance/popGrowthParamsMCMC.rds")


####
#### Organize parameter values ---------------------------------
####
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=4))
pGrow2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("betaSpp", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intercept", times=4*3000),
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
  growNowYr <- subset(pGrowYrs, Year==doYear)
  growNowYr <- subset(growNowYr, Iter==doNow)
  growNowYr <- subset(growNowYr, Spp==sppSim)
  iID <- which(growNowYr$Coef=="intYr")
  intercept <- growNowYr$value[iID]
  gID <- which(growNow$Coef=="gInt")
  intG <- growNow$value[gID]
  sID <- which(growNowYr$Coef=="beta")
  size <- growNowYr$value[sID]
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  newN <- intercept+size*N+sum(climEffs*climate)
  newN <- rpois(1,exp(newN))
  return(newN)
}


####
#### Run simulations -----------------------------------------------------
####
outD <- data.frame(year=NA, cover=NA, sim=NA, species=NA)

for(i in 1:length(sppList)){
  sppSim <- sppList[i]
  nSim <- 2
  yearsN <- length(unique(allD$year))
  years <- unique(allD$year)+1900
  yearsID <- unique(allD$year)
  Nsave <- matrix(ncol=yearsN, nrow=nSim)
  sppD <- subset(allD, Species==sppSim)
  Nsave[,1] <- round(mean(subset(sppD, year==yearsID[1])$percCover)*10000)
  
  for(sim in 1:nSim){
    for(yr in 2:yearsN){
      N <- Nsave[sim,yr-1]
      climate <- subset(climD, year==years[yr-1])[,c(3,5,4,6)]
      Nout <- growFunc(pGrow=pGrowAll, pGrowYrs=pGrowYrs, N=N, climate=climate, simsPerYear=length(NforG), doYear=years[yr], sppSim=sppSim)
      Nsave[sim,yr] <- Nout
      print(c(sim, yr, sppSim))
    }#end year loop
  }#end sim loop
  
  dN <- as.data.frame(Nsave)
  colnames(dN) <- years
  nM <- melt(dN)
  nM$sim <- rep(1:nSim, length(years))
  nM$species <- rep(sppSim, nSim*length(years))
  colnames(nM)[1:2] <- c("year", "cover")
  outD <- rbind(outD, nM)
}

####
#### Make output data frame ------------------------------------------------
####
# noNA <- which(is.na(outD$cover)!=TRUE)
outD <- outD[2:nrow(outD),]
quadD <- ddply(allD, .variables = c("year", "Species"), .fun = summarise,
               year = mean(year),
               cover = mean(percCover))
quadD$year <- quadD$year+1900
colnames(quadD) <- c("species", "year", "obsCover")
plotD <- merge(outD, quadD)
d1 <- subset(plotD, species==sppList[1])
d2 <- subset(plotD, species==sppList[2])
d3 <- subset(plotD, species==sppList[3])
d4 <- subset(plotD, species==sppList[4])

g1 <- ggplot(data=d1)+
#   geom_line(aes(x=year, y=cover/100, group=sim), alpha=1, color="purple")+
  geom_line(aes(x=year, y=obsCover*100, group=NA), color="grey25")+
  geom_point(aes(x=year, y=obsCover*100), size=4, color="grey25")+
  ylab("Mean Cover (%)")+
  xlab("Year")+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle(sppList[1])
g2 <- g1 %+% d2 + ggtitle(sppList[2])
g3 <- g1 %+% d3 + ggtitle(sppList[3])
g4 <- g1 %+% d4 + ggtitle(sppList[4])

g <- arrangeGrob(g1,g2,g3,g4)
g
