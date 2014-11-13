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
                    colonizes= NA,
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
  
  #Colonization observations
  colD <- subset(sppD,lag.cover==0)
  colD$colonizes <- ifelse(colD$totCover>0,1,0)
  colD$yearID <- colD$year #for random year offset on intercept
  colD$group <- substring(colD$quad, 1, 1)
  backD <- rbind(backD, colD)
}#end species loop

colD <- backD[2:nrow(backD),]

####
#### Set up data structure for JAGS
####
nGrp <- length(unique(colD$group))
nYrs <- length(unique(colD$year))
nObs <- nrow(colD)
C <- colD$colonizes
yrs <- (colD$year)-32
grp <- as.numeric(as.factor(colD$group))
TmeanSpr1 <- colD$TmeanSpr1
TmeanSpr2 <- colD$TmeanSpr2
ppt1 <- colD$ppt1
ppt2 <- colD$ppt2
nSpp <- length(unique(colD$Species))
spp <- as.numeric(as.factor(colD$Species))

####
#### Run MCMC
####
iterations <- 50000
adapt <- 10000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, yrs=yrs, grp=grp,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, spp=spp, nSpp=nSpp)
mod <- jags.model("colonizationAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("intercept", "intYr", "intG", "temp1", "temp2", "rain1", "rain2"),
                    n.iter=iterations, n.thin=10)
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

####
#### Check for convergence
####
gelmDiag <- gelman.diag(out)
# heidel.diag(out)
# gelman.plot(out)

pdf("colonizationOutPlots.pdf")
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

sppNames <- c(rep("all", nGrp), rep(sppList, 5))
outStat$species <- sppNames
outQuant$species <- sppNames

saveRDS(outC, file = "colonizationParamsMCMC.rds")
write.csv(gelmDiag[[1]], file="colonizationGelman.csv")
write.csv(outStat, file="colonizationStats.csv")
write.csv(outQuant, file="colonizationQuants.csv")
write.csv(outDeviance, file="colonizationDeviance.csv")




