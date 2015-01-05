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
gPars <- as.data.frame(read.csv("growthStats.csv"))
years <- unique(growD$year)
gPars$year <- c(rep(years, each=4),
                rep(NA, times=4*6),
                rep(years, each=4),
                rep(NA, times=4*5))
gPars$Group <- c(rep(NA, times=13*4),
                 rep(sort(unique(growD$Group)), each=4),
                 rep(NA, times=4*13),
                 rep(NA, times=4*5))

####
#### OK, do predictions
####
y <- numeric(nrow(growD))
for(i in 1:nrow(growD)){
  gNow <- as.numeric(growD$Group[i])
  sppNow <- growD$species[i]
  yrNow <- growD$year[i]
  beta <- subset(gPars[grep("beta", x = gPars$X),], species==sppNow & year==yrNow)$Mean
  int <- subset(gPars[grep("intYr", x = gPars$X),], species==sppNow & year==yrNow)$Mean
  intG <- subset(gPars[grep("intG", x = gPars$X),], species==sppNow & Group==gNow)$Mean
  nb <- subset(gPars[grep("nb", x = gPars$X),], species==sppNow)$Mean
  rain1 <- subset(gPars[grep("rain1", x = gPars$X),], species==sppNow)$Mean
  rain2 <- subset(gPars[grep("rain2", x = gPars$X),], species==sppNow)$Mean
  temp1 <- subset(gPars[grep("temp1", x = gPars$X),], species==sppNow)$Mean
  temp2 <- subset(gPars[grep("temp2", x = gPars$X),], species==sppNow)$Mean
  x <- log(growD$area.t0[i])
  crowd <- growD$crowd[i]
  TmeanSpr1 <- growD$TmeanSpr1[i]
  TmeanSpr2 <- growD$TmeanSpr2[i]
  ppt1 <- growD$ppt1[i]
  ppt2 <- growD$ppt2[i]
  y[i] <- int+intG+beta*x+nb*crowd+rain1*ppt1+rain2*ppt2+temp1*TmeanSpr1+temp2*TmeanSpr2
  print(paste(i, "of", nrow(growD)))
}

yR2 <- (y-log(growD$area.t0))^2
# plot(y, yR2)
growD$yHat <- y
growD$yR2 <- yR2
varPars <- matrix(nrow=length(sppList), ncol=2)
for(i in 1:length(sppList)){ 
  outVar <- nls(yR2~a*exp(b*yHat),start=list(a=1,b=0), data=subset(growD, species==sppList[i]))
  varPars[i,] <- coef(outVar)
}

varPars <- as.data.frame(varPars)
colnames(varPars) <- c("a", "b")
varPars$species <- sppList

####
#### Write results to file
####
saveRDS(varPars, "varianceParams.rds")
  
  