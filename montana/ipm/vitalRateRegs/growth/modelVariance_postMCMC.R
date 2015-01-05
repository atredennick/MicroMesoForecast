#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#### Get variance parameters for growth regression (IPM) ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#Makes variance a function of plant size
#example from Peter's Idaho code: 
##### fit variance                                #####
##### x=fitted(out)                                #####
##### y=resid(out)^2                               #####
##### plot(x,y)                                    #####
##### outVar=nls(y~a*exp(b*x),start=list(a=1,b=0)) #####

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(plyr)
library(reshape2)

sppList=sort(c("BOGR","HECO","PASM","POSE"))

####
#### Read in data by species and make one long data frame
####
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area.t1=NA,
                   area.t0=NA,
                   age=NA,
                   allEdge=NA,
                   distEdgeMin=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
  sppD$species <- doSpp
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]
climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

crowd <- c(read.csv("BOGRgrowthCrowding.csv")[,2], 
           read.csv("HECOgrowthCrowding.csv")[,2],
           read.csv("PASMgrowthCrowding.csv")[,2],
           read.csv("POSEgrowthCrowding.csv")[,2])
growD$crowd <- crowd

####
#### Bring in mean MCMC-based parameters for predictions
####
gPars <- read.csv("growthStats.csv")
gPars$Year <- 
  
  
  