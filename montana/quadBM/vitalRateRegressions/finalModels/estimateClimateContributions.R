#This script fits climate and constant models to estimate contribution of climate variables to the model

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(INLA)

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")

files <- list.files()
files <- files[grep("growth", files)]

random <- "+f(yearID, model='iid', prior='normal',param=c(0,0.001))+
              f(year, percLagCover, model='iid', prior='normal', param=c(0,0.001))+
              f(group, model='iid', prior='normal',param=c(0,0.001))"
randomNoTemp <- "+f(yearID, model='iid', prior='normal',param=c(0,0.001))+
                  f(group, model='iid', prior='normal',param=c(0,0.001))"

SSRs <- matrix(ncol=5, nrow=(length(sppList)))

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  tmpFile <- files[grep(doSpp, files)]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  
  #load full model
  fullMod <- readRDS(tmpFile)
  
  #fit constant and climate models
  constant <- as.formula(paste("percCover ~ percLagCover", randomNoTemp)) 
  constantMod <- inla(constant, data=growD,
                       family=c("beta"), verbose=TRUE,
                       control.compute=list(dic=T,mlik=T),
                       control.predictor = list(link = 1),
                       control.inla = list(h = 1e-10))
  fixedCovs <- fullMod$names.fixed[2:length(fullMod$names.fixed)]
  fixedCovs <- paste(fixedCovs, collapse="+")
  tmpFix <- paste("percCover ~ ", fixedCovs, sep="")
  climate <- as.formula(paste(tmpFix, randomNoTemp)) 
  climateMod <- inla(climate, data=growD,
                       family=c("beta"), verbose=TRUE,
                       control.compute=list(dic=T,mlik=T),
                       control.predictor = list(link = 1),
                       control.inla = list(h = 1e-10))
  
  ssrFull <- fullMod$dic$mean.deviance
  ssrConstant <- constantMod$dic$mean.deviance
  ssrClimate <- climateMod$dic$mean.deviance
  
  climContribution <- (ssrClimate-ssrConstant)/(ssrFull-ssrConstant)
  SSRs[spp,] <- c(length(growD$percCover), ssrConstant, ssrClimate, ssrFull, climContribution) 

}#end species loop
