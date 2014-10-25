#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(reshape2)
library(plyr)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")

doSpp <- sppList[1]

#load vital rate parameters
files <- list.files("../vitalRateRegressions/finalModels/vitalRateParameters/")
paramFiles <- files[grep(doSpp, files)]
pCol <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[1], sep=""))
pGrow <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[2], sep=""))
pSurv <- readRDS(paste("../vitalRateRegressions/finalModels/vitalRateParameters/", paramFiles[3], sep=""))


####
#### Vital rate functions -----------------------------------------------
####
survFunc <- function(pSurv, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pSurv$fixed[1,1], pSurv$fixed[1,2])
  size <- rnorm(simsPerYear, pSurv$fixed[2,1], pSurv$fixed[2,2])
  yearID <- which(pSurv$yearInt$ID==doYear)
  yearOff <- rnorm(simsPerYear, pSurv$yearInt[yearID,2], pSurv$yearInt[yearID,3])
  yearSize <- rnorm(simsPerYear, pSurv$yearSize[yearID,2], pSurv$yearSize[yearID,3])
  
  newN <- intercept+yearOFF+(size+yearSize)*N
  newN <- 1/(1+newN)
  return(newN)
}

growFunc <- function(pGrow, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pGrow$fixed[1,1], pGrow$fixed[1,2])
  size <- rnorm(simsPerYear, pGrow$fixed[2,1], pGrow$fixed[2,2])
  yearID <- which(pGrow$yearInt$ID==doYear)
  yearOff <- rnorm(simsPerYear, pGrow$yearInt[yearID,2], pGrow$yearInt[yearID,3])
  yearSize <- rnorm(simsPerYear, pGrow$yearSize[yearID,2], pGrow$yearSize[yearID,3])
  
  newN <- intercept+yearOFF+(size+yearSize)*N
  newN <- 1/(1+newN)
  return(newN)
}

colFunc <- function(pCol, N, climate, simsPerYear, doYear){
  intercept <- rnorm(simsPerYear, pCol$fixed[1,1], pCol$fixed[1,2])
  size <- rnorm(simsPerYear, pCol$fixed[2,1], pCol$fixed[2,2])
  yearID <- which(pCol$yearInt$ID==doYear)
  yearOff <- rnorm(simsPerYear, pCol$yearInt[yearID,2], pCol$yearInt[yearID,3])
  yearSize <- rnorm(simsPerYear, pCol$yearSize[yearID,2], pCol$yearSize[yearID,3])
  
  newN <- intercept+yearOFF+(size+yearSize)*N
  newN <- 1/(1+newN)
  return(newN*0.01)
}


####
#### Run simulations -----------------------------------------------------
####






