#This script fits climate and constant models to estimate contribution of climate variables to the model

####
#### FOR SURVIVAL REGRESSIONS
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
  
  #Survival observations
  survD <- subset(sppD,lag.cover>0)
  survD$survives <- ifelse(survD$totCover>0,1,0)
  survD$yearID <- survD$year #for random year offset on intercept
  survD$group <- substring(survD$quad, 1, 1)
  survD$percCover <- survD$totCover/10000
  survD$percLagCover <- survD$lag.cover/10000
  
  #load full model
  fullMod <- readRDS(tmpFile)
  
  #fit constant and climate models
  constant <- as.formula(paste("survives ~ percLagCover", randomNoTemp)) 
  constantMod <- inla(constant, data=survD,
                      family=c("binomial"), verbose=TRUE,
                      control.compute=list(dic=T,mlik=T),
                      control.inla = list(h = 1e-10),
                      control.predictor = list(link = 1),
                      Ntrials=rep(1,nrow(survD)))
  
  fixedCovs <- fullMod$names.fixed[2:length(fullMod$names.fixed)]
  fixedCovs <- paste(fixedCovs, collapse="+")
  tmpFix <- paste("survives ~ ", fixedCovs, sep="")
  climate <- as.formula(paste(tmpFix, randomNoTemp)) 
  climateMod <- inla(climate, data=survD,
                     family=c("binomial"), verbose=TRUE,
                     control.compute=list(dic=T,mlik=T),
                     control.inla = list(h = 1e-10),
                     control.predictor = list(link = 1),
                     Ntrials=rep(1,nrow(survD)))
  
  ssrFull <- fullMod$dic$mean.deviance
  ssrConstant <- constantMod$dic$mean.deviance
  ssrClimate <- climateMod$dic$mean.deviance
  
  climContribution <- (ssrClimate-ssrConstant)/(ssrFull-ssrConstant)
  SSRs[spp,] <- c(length(survD$percCover), ssrConstant, ssrClimate, ssrFull, climContribution) 

}#end species loop

ssrTab <- data.frame(Species = sppList,
                     N = SSRs[,1],
                     SSR1 = SSRs[,2],
                     SSR2 = SSRs[,3],
                     SSR3 = SSRs[,4],
                     SSRc = SSRs[,5])
colnames(ssrTab) <- c("Species", "N", "Constant model", "Climate model", 
                      "Full model", "Contribution of climate covariates")
write.csv(ssrTab, "survivalClimateContributionTable.csv", row.names=FALSE)

