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
allD <- read.csv("../../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../../weather/Climate.csv")
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
TmeanSpr2 <- growD$TmeanSpr2
ppt1 <- growD$ppt1
ppt2 <- growD$ppt2
nSpp <- length(unique(growD$Species))
spp <- as.numeric(as.factor(growD$Species))
X <- growD$percLagCover

####
#### Run MCMC
####
iterations <- 50000
adapt <- 10000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, X=X, yrs=yrs, grp=grp,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, spp=spp, nSpp=nSpp)
mod <- jags.model("growthAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("intYr", "intercept", "beta", "betaSpp", "intG", "temp1", "temp2", "rain1", "rain2", "tau"),
                    n.iter=iterations, n.thin=10)
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

####
#### Check for convergence
####
gelmDiag <- gelman.diag(out)
# heidel.diag(out)
# gelman.plot(out)

pdf("growthOutPlots.pdf")
plot(out, auto.layout=FALSE)
dev.off()

####
#### Convert to dataframe for export and get other summaries
####
outC <- rbind(out[[1]][(iterations-999):iterations,], 
              out[[2]][(iterations-999):iterations,], 
              out[[3]][(iterations-999):iterations,])

outStat <- as.data.frame(summary(out)$stat)
outQuant <- as.data.frame(summary(out)$quantile)
outDeviance <- as.data.frame(summary(dic$deviance, mean)$stat)

sppNames <- c(rep(sppList, 13+6+13+4+1+1+1))
outStat$species <- sppNames
outQuant$species <- sppNames

saveRDS(outC, file = "growthParamsMCMC.rds")
write.csv(gelmDiag[[1]], file="growthGelman.csv")
write.csv(outStat, file="growthStats.csv")
write.csv(outQuant, file="growthQuants.csv")
write.csv(outDeviance, file="growthDeviance.csv")
# 
# 
# 
# 
