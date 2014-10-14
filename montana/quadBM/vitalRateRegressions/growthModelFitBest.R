#This script pulls in the DICs and then selects/fits the best model for each species

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#bring in DIC values
dicD <- read.csv("growthDIC.csv")

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

#go species by species, sort DICs, save best model
deltaD <- data.frame(modelNum = NA, Species=NA, deltaDIC=NA)
for(i in 1:length(sppList)){
  dicTmp <- subset(dicD, Species==sppList[i])
  dicTmp <- dicTmp[with(dicTmp, order(DIC)), ]
  dicTmp$deltaDIC <- NA
  for(j in 1:nrow(dicTmp)){
    dicTmp$deltaDIC[j] <- dicTmp$DIC[j] - dicTmp$DIC[1]
  }
  subTmp <- subset(dicTmp, deltaDIC<5)
  subTmp <- subTmp[,c(2,3,5)]
  deltaD <- rbind(deltaD, subTmp)
}

#write the top models for each species
write.csv(deltaD[2:nrow(deltaD),], "topModels_Growth.csv", row.names=FALSE)

