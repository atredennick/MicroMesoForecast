#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)
library(parallel) 
library(snowfall) 
library(rlecuyer) 
load.module("dic")
source("mont_qbm_jags.R")

cps=detectCores()
sfInit(parallel=TRUE, cpus=cps)
sfExportAll()
sfClusterSetupRNG()

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

backD <- data.frame(climYear=NA,
                    quad = NA,
                    year= NA,
                    totCover= NA,
                    Species= NA,
                    propCover= NA,
                    lag.cover= NA,
                    pptLag= NA,
                    ppt1= NA,
                    TmeanSpr1= NA,
                    ppt2= NA,
                    TmeanSpr2= NA,
                    TmeanSum1= NA,
                    TmeanSum2= NA,
                    yearID= NA,
                    group = NA,
                    percCover = NA,
                    percLagCover = NA)

#loop through species and remake data frame
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
  sppD<-merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  backD <- rbind(backD, growD)
}#end species loop
growD <- backD[2:nrow(backD),]
 

####
#### Setup climate covariates
####
growD$interaction1 <- with(growD, ppt1*TmeanSpr1)
growD$interaction2 <- with(growD, ppt2*TmeanSpr2)

grow_now <- subset(growD, Species=="BOGR")
X <- list(X1 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction1", "interaction2")],
          X2 = grow_now[,c("ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction1", "interaction2")],
          X3 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction1")],
          X4 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction2")],
          X5 = grow_now[,c("ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction1")],
          X6 = grow_now[,c("ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction2")],
          X7 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")],
          X8 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "interaction1")],
          X9 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr2", "interaction2")],
          X10 = grow_now[,c("pptLag", "ppt1", "TmeanSpr1", "TmeanSpr2", "interaction1")],
          X11 = grow_now[,c("pptLag", "ppt2", "TmeanSpr1", "TmeanSpr2", "interaction2")],
          X12 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1")],
          X13 = grow_now[,c("pptLag", "ppt1", "ppt2", "TmeanSpr2")],
          X14 = grow_now[,c("pptLag", "ppt1", "TmeanSpr1", "TmeanSpr2")],
          X15 = grow_now[,c("pptLag", "ppt2", "TmeanSpr1", "TmeanSpr2")],
          X16 = grow_now[,c("ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")])

####
#### Send to model fitting function
####
iters <- 50000
dic <- numeric(length(X))
for(i in 1:length(X)){
  model_dic <- mont_qbm_jags(Y = grow_now$percCover, size = grow_now$percLagCover, X = X[[i]], 
                             groups = as.numeric(as.factor(grow_now$group)), 
                             species = as.numeric(as.factor(grow_now$Species)),
                             years = (grow_now$year-32), iters=iters, adapt.iters=500, thins=5)
  dic[i] <- sum(model_dic$deviance)+sum(model_dic$penalty)
}


