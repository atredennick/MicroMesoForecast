#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)
load.module("dic")

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
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

####
#### Set up data for JAGS model
####
dataJ <- list(Y = log(growD$area.t1),
              X = log(growD$area.t0),
              nObs = nrow(growD),
              grp = as.numeric(growD$Group),
              nGrp = length(unique(growD$Group)),
              spp = as.numeric(as.factor(growD$species)),
              nSpp = length(unique(growD$species)),
              yrs = (growD$year - 31),
              nYrs = length(unique(growD$year)),
              TmeanSpr1 = growD$TmeanSpr1,
              TmeanSpr2 = growD$TmeanSpr2,
              ppt1 = growD$ppt1,
              ppt2 = growD$ppt2)

####
#### Run MCMC from JAGS
####
iterations <- 5000
adapt <- 1000
mod <- jags.model("growthAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("intYr", "beta", "intG", "temp1", "temp2", "rain1", "rain2"),
                    n.iter=iterations, n.thin=10)
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

plot(out)


