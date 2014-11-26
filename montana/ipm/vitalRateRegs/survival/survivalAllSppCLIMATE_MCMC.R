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
                   area=NA,
                   survives=NA,
                   age=NA,
                   distEdgeMin=NA,
                   allEdge=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- read.csv(paste("../../../speciesData/", doSpp, "/survD.csv", sep=""))
  sppD$species <- doSpp
  outD <- rbind(outD, sppD)
}

survD <- outD[2:nrow(outD),]
climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
survD <- merge(survD,climD)
survD$Group=as.factor(substr(survD$quad,1,1))

####
#### Set up data for JAGS model
####
dataJ <- list(Y = survD$survives,
              X = log(survD$area),
              nObs = nrow(survD),
              grp = as.numeric(survD$Group),
              nGrp = length(unique(survD$Group)),
              spp = as.numeric(as.factor(survD$species)),
              nSpp = length(unique(survD$species)),
              yrs = (survD$year - 31),
              nYrs = length(unique(survD$year)),
              crowd = crowd,
              TmeanSpr1 = survD$TmeanSpr1,
              TmeanSpr2 = survD$TmeanSpr2,
              ppt1 = survD$ppt1,
              ppt2 = survD$ppt2)

####
#### Run MCMC from JAGS
####
iterations <- 50000
adapt <- 10000
mod <- jags.model("survivalAllSppCLIMATE_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

outDeviance <- as.data.frame(summary(dic$deviance, mean)$stat)
write.csv(outDeviance, file="survivalDevianceCLIMATE.csv")


