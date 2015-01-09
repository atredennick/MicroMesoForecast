#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#### Script to fit population growth model to aggregated Montana individual-level data ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

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
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../weather/Climate.csv")
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
                    group = NA)

#loop through species and remake data frame
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
  tmp=sppD[,c("quad","year","propCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD2=merge(sppD,climD,by.x="climYear",by.y="year")
  sppD2$yearID <- sppD2$year #for random year offset on intercept
  sppD2$group <- substring(sppD2$quad, 1, 1)
  backD <- rbind(backD, sppD2)
}#end species loop

modD <- backD[2:nrow(backD),]
naID <- which(is.na(modD$lag.cover)==FALSE)
modD <- modD[naID,]
####
#### Set up data structure for JAGS
####
nGrp <- length(unique(modD$group))
nYrs <- length(unique(modD$year))
nObs <- nrow(modD)
C <- round(modD$propCover*100)
yrs <- (modD$year)-32
grp <- as.numeric(as.factor(modD$group))
TmeanSpr1 <- modD$TmeanSpr1
TmeanSpr2 <- modD$TmeanSpr2
ppt1 <- modD$ppt1
ppt2 <- modD$ppt2
nSpp <- length(unique(modD$Species))
spp <- as.numeric(as.factor(modD$Species))
X <- round(modD$lag.cover*100)

####
#### Run MCMC
####
iterations <- 50000
adapt <- 10000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, X=X, yrs=yrs, grp=grp,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, spp=spp, nSpp=nSpp)
mod <- jags.model("popAbundanceModelJAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("intYr", "beta", "intG", "temp1", "temp2", "rain1", "rain2"),
                    n.iter=iterations, n.thin=10)
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

####
#### Check for convergence
####
gelmDiag <- gelman.diag(out)
# heidel.diag(out)
# gelman.plot(out)
# 
pdf("survivalOutPlots.pdf")
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

sppNames <- rep(sppList, 1+5+6)
outStat$species <- sppNames
outQuant$species <- sppNames

saveRDS(outC, file = "popGrowthParamsMCMC.rds")
write.csv(gelmDiag[[1]], file="popGrowthGelman.csv")
write.csv(outStat, file="popGrowthStats.csv")
write.csv(outQuant, file="popGrowthQuants.csv")
write.csv(outDeviance, file="popGrowthDeviance.csv")
