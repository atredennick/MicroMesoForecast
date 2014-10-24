#This script pulls in the DICs and then selects/fits the best model for each species

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(INLA)

#bring in DIC values
dicD <- read.csv("topModels_Growth.csv")

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../weather/Climate.csv")

#get all the models and get best from DICs for each species
source("fixedEffectsModels.R")


####
#### Get best models for each species (by number ID)
####

bestMod <- data.frame(Species = sppList,
                      BestMod = rep(NA, length(sppList)))
for(i in 1:length(sppList)){
  tmpDIC <- subset(dicD, Species==sppList[i])
  bestMod[i,2] <- tmpDIC[1,1]
}

####
#### Fit best model for each species
####

#set up unchanging random effects
random <- "+f(yearID, model='iid', prior='normal',param=c(0,0.001))+
              f(year, percLagCover, model='iid', prior='normal', param=c(0,0.001))+
              f(group, model='iid', prior='normal',param=c(0,0.001))"

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
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
  
  #fit model
  modNum <- bestMod[spp,2]
  fixed <- fixed.forms[[modNum]]
  formula <- as.formula(paste(fixed, random)) 
  out <- inla(formula, data=growD,
               family=c("beta"), verbose=TRUE,
               control.compute=list(dic=T,mlik=T),
               control.inla = list(h = 1e-10))
  saveRDS(out, file = paste("./finalModels/growth", doSpp, ".rds", sep=""))
}#end species loop






