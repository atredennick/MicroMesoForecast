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
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

crowd <- c(read.csv("BOGRgrowthCrowding.csv")[,2], 
           read.csv("HECOgrowthCrowding.csv")[,2],
           read.csv("PASMgrowthCrowding.csv")[,2],
           read.csv("POSEgrowthCrowding.csv")[,2])


# library(lme4)
# D <- subset(growD, species=="BOGR")
# D$logarea.t1 <- log(D$area.t1)
# D$logarea.t0 <- log(D$area.t0)
# outlm=lmer(logarea.t1~logarea.t0+ppt1+TmeanSpr1+ppt2+TmeanSpr2+
#            (1|Group)+(logarea.t0|year),data=D)

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
              crowd = crowd,
              TmeanSpr1 = growD$TmeanSpr1,
              TmeanSpr2 = growD$TmeanSpr2,
              ppt1 = growD$ppt1,
              ppt2 = growD$ppt2)

####
#### Run MCMC from JAGS
####
iterations <- 50000
adapt <- 10000
mod <- jags.model("growthAllSppCONSTANT_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

####
#### Convert to dataframe for export and get other summaries
####
outDeviance <- as.data.frame(summary(dic$deviance, mean)$stat)
write.csv(outDeviance, file="growthDevianceCONSTANT.csv")








