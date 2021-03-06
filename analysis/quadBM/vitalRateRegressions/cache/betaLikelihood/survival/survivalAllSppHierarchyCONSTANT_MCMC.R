#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)
load.module("dic")

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)

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
                    survives= NA,
                    yearID= NA,
                    group = NA)

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
  
  #Survival observations
  survD <- subset(sppD,lag.cover>0)
  survD$survives <- ifelse(survD$totCover>0,1,0)
  survD$yearID <- survD$year #for random year offset on intercept
  survD$group <- substring(survD$quad, 1, 1)
  backD <- rbind(backD, survD)
}#end species loop

survD <- backD[2:nrow(backD),]

####
#### Set up data structure for JAGS
####
nGrp <- length(unique(survD$group))
nYrs <- length(unique(survD$year))
nObs <- nrow(survD)
C <- survD$survives
yrs <- (survD$year)-32
grp <- as.numeric(as.factor(survD$group))
TmeanSpr1 <- survD$TmeanSpr1
TmeanSpr2 <- survD$TmeanSpr2
ppt1 <- survD$ppt1
ppt2 <- survD$ppt2
nSpp <- length(unique(survD$Species))
spp <- as.numeric(as.factor(survD$Species))
X <- log(survD$lag.cover)

####
#### Run MCMC
####
iterations <- 50000
adapt <- 10000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, X=X, yrs=yrs, grp=grp,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, spp=spp, nSpp=nSpp)
mod <- jags.model("survivalAllSppCONSTANT_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

####
#### Convert to dataframe for export and get other summaries
####
outDeviance <- as.data.frame(summary(dic$deviance, mean)$stat)

write.csv(outDeviance, file="survivalDevianceCONSTANT.csv")




