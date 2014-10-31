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

#loop through species and fit models
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
  colD$percCover <- colD$totCover/10000
  colD$percLagCover <- colD$lag.cover/10000
  
  
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
  
  ####
  #### Run MCMC
  ####
  iterations <- 50000
  adapt <- 30000
  dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, yrs=yrs, grp=grp,
                TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2)
  mod <- jags.model("colonizationAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
  update(mod, n.iter = (iterations/2))
  out <- coda.samples(mod, c("interceptMu"),
                      n.iter=iterations, n.thin=10)
  dic <- jags.samples(mod, c("deviance"),
                      n.iter=iterations, n.thin=10)
  
  ####
  #### Check for convergence
  ####
  gelmDiag <- gelman.diag(out)
  # heidel.diag(out)
  # gelman.plot(out)
  
  pdf(paste(doSpp, "_growthOutPlots.pdf", sep=""))
  plot(out, auto.layout=FALSE)
  dev.off()
  
  ####
  #### Convert to dataframe for export and get other summaries
  ####
  outC <- rbind(out[[1]][(iterations-999):iterations,], 
                out[[2]][(iterations-999):iterations,], 
                out[[3]][(iterations-999):iterations,])
  
  outStat <- summary(out)$stat
  outQuant <- summary(out)$quantile
  outDeviance <- summary(dic$deviance, mean)$stat
  
  saveRDS(outC, file = paste(doSpp, "_SurvivalParamsMCMC.rds", sep=""))
  write.csv(gelmDiag[[1]], file=paste(doSpp, "_survivalGelman.csv", sep=""))
  write.csv(outStat, file=paste(doSpp, "_survivalStats.csv", sep=""))
  write.csv(outQuant, file=paste(doSpp, "_survivalQuants.csv", sep=""))
  write.csv(outDeviance, file=paste(doSpp, "_survivalDeviance.csv", sep=""))
}





