#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)
library(ggplot2)
library(plyr)
library(reshape2)
library(ggthemes)
load.module("dic")

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale=TRUE)

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
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
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
#### Set up data structure for JAGS
####
nGrp <- length(unique(growD$group))
nYrs <- length(unique(growD$year))
nObs <- nrow(growD)
C <- growD$percCover
yrs <- (growD$year)-32
grp <- as.numeric(as.factor(growD$group))
TmeanSpr1 <- growD$TmeanSpr1
TmeanSpr1V <- growD$TmeanSpr1
TmeanSpr2 <- growD$TmeanSpr2
ppt1 <- growD$ppt1
ppt2 <- growD$ppt2
nSpp <- length(unique(growD$Species))
spp <- as.numeric(as.factor(growD$Species))
sppId1 <- which(spp==1)
sppId2 <- which(spp==2)
sppId3 <- which(spp==3)
sppId4 <- which(spp==4)
X <- growD$percLagCover

####
#### Run MCMC
####
iterations <- 5000
adapt <- 100
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, X=X, yrs=yrs, grp=grp, TmeanSpr1V=TmeanSpr1V,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, spp=spp, nSpp=nSpp)
mod <- jags.model("growthAllSpp_JAGS_ClimateTrack.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("ppt1P"),
                    n.iter=iterations, n.thin=10)

outQuant <- as.data.frame(summary(out)$quantile)
outQuant$spp <- spp 
outQuant$yr <- yrs+1932-1 
colnames(outQuant)[1:5] <- c("lo2", "lo25", "med", "up75", "up97")
outD <- ddply(outQuant, .(yr, spp), summarise,
              avgT = mean(med),
              lowT = mean(lo2),
              highT = mean(up97))
outD <- merge(outD, climD, by.x = "yr", by.y = "year")
ggplot(outD, aes(x=yr, y=avgT, group=as.factor(spp),  color=as.factor(spp)))+
#   geom_line(size=1)+
  geom_line(aes(x=yr, y=lowT), linetype=2, size=1)+
  geom_line(aes(x=yr, y=highT), linetype=2, size=1)+
  geom_line(aes(x=yr, y=ppt1), color="grey25", size=1)+
  facet_grid(spp~.)+
  guides(color=FALSE)+
  scale_color_manual(values=c("darkorange", "steelblue", "darkgreen", "purple"))+
  theme_classic()