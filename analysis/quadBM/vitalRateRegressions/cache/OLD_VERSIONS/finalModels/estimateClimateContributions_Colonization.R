#This script fits climate and constant models to estimate contribution of climate variables to the model

####
#### FOR COLONIZATION REGRESSIONS
####

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
files <- files[grep("survival", files)]

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
  
  #Colonization observations
  colD <- subset(sppD,lag.cover==0)
  colD$colonizes <- ifelse(colD$totCover>0,1,0)
  colD$yearID <- colD$year #for random year offset on intercept
  colD$group <- substring(colD$quad, 1, 1)
  colD$percCover <- colD$totCover/10000
  colD$percLagCover <- colD$lag.cover/10000
  
  #load full model
  fullMod <- readRDS(tmpFile)
  
  #fit constant and climate models
  constant <- as.formula(paste("colonizes ~ ", randomNoTemp)) 
  constantMod <- inla(constant, data=colD,
                      family=c("binomial"), verbose=TRUE,
                      control.compute=list(dic=T,mlik=T),
                      control.inla = list(h = 1e-10),
                      control.predictor = list(link = 1),
                      Ntrials=rep(1,nrow(colD)))
  
  fixedCovs <- fullMod$names.fixed[2:length(fullMod$names.fixed)]
  fixedCovs <- paste(fixedCovs, collapse="+")
  tmpFix <- paste("colonizes ~ ", fixedCovs, sep="")
  climate <- as.formula(paste(tmpFix, randomNoTemp)) 
  climateMod <- inla(climate, data=colD,
                     family=c("binomial"), verbose=TRUE,
                     control.compute=list(dic=T,mlik=T),
                     control.inla = list(h = 1e-10),
                     control.predictor = list(link = 1),
                     Ntrials=rep(1,nrow(colD)))
  
  ssrFull <- fullMod$dic$mean.deviance
  ssrConstant <- constantMod$dic$mean.deviance
  ssrClimate <- climateMod$dic$mean.deviance
  
  climContribution <- (ssrClimate-ssrConstant)/(ssrFull-ssrConstant)
  SSRs[spp,] <- c(length(colD$percCover), ssrConstant, ssrClimate, ssrFull, climContribution) 

}#end species loop

ssrTab <- data.frame(Species = sppList,
                     N = SSRs[,1],
                     SSR1 = SSRs[,2],
                     SSR2 = SSRs[,3],
                     SSR3 = SSRs[,4],
                     SSRc = SSRs[,5])
colnames(ssrTab) <- c("Species", "N", "Constant model", "Climate model", 
                      "Full model", "Contribution of climate covariates")
write.csv(ssrTab, "colonizationClimateContributionTable.csv", row.names=FALSE)

