####
#### Check climate covariates and random year effect correlations
####

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(plyr)
library(reshape2)

climD <- read.csv("../../weather/Climate.csv")

####
#### Growth
####
gParms <- read.csv("./growth/growthStats.csv")
bID <- grep("beta", x = gParms$X)
iID <- grep("intYr", x = gParms$X)
gSub <- gParms[c(bID, iID),]
gSub$coef <- c(rep("beta", 52), rep("intercept", 52))
gSub$year <- rep(rep(climD$year, each=4), times=2)
gSub <- merge(gSub, climD[,c(1,3,4,5,6)], by="year")
gBeta <- subset(gSub, coef=="beta")
gInt <- subset(gSub, coef=="intercept")
cor(gBeta[,3], gBeta[,9:12])
cor(gInt[,3], gInt[,9:12])
