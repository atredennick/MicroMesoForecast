#Script to get parameters for final growth models for each species and export for IBM simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(INLA)

####
#### Load final models by species and vital rate
####

vitals <- c("growth", "survival", "colonization")

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

files <- list.files()

for(i in 1:length(vitals)){
  fileTmp <- files[grep(glob2rx(paste(vitals[i], "*.rds", sep="")), files)]
  for(j in 1:length(sppList)){
    fileTmp2 <- fileTmp[grep(sppList[j], fileTmp)]
    modTmp <- readRDS(fileTmp2)
    
    fixed <- modTmp$summary.fixed[,1:2]
    yearInt <- modTmp$summary.random$yearID[,1:3]
    yearSize <- modTmp$summary.random$year[,1:3]
    groupInt <- modTmp$summary.random$group[,1:3]
    
    params <- list(fixed=fixed, yearInt=yearInt, yearSize=yearSize, groupInt=groupInt)
    saveRDS(params, file = paste("./vitalRateParameters/", vitals[i], sppList[j], "_params.rds", sep=""))
    
  }#end species loop
}#end vital rate loop

