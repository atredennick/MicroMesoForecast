#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
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
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  
  
  ####
  #### Set up data structure for JAGS
  ####
  nGrp <- length(unique(growD$group))
  nYrs <- length(unique(growD$year))
  nObs <- nrow(growD)
  size <- growD$percLagCover
  Y <- growD$percCover
  yrs <- (growD$year)-32
  grp <- as.numeric(as.factor(growD$group))
  TmeanSpr1 <- growD$TmeanSpr1
  TmeanSpr2 <- growD$TmeanSpr2
  ppt1 <- growD$ppt1
  ppt2 <- growD$ppt2
  
  ####
  #### Run MCMC
  ####
  iterations <- 2000
  dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, size=size, Y=Y, yrs=yrs, grp=grp,
                TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2)
  mod <- jags.model(paste("growth", doSpp, "_JAGS.R", sep=""), data=dataJ, n.chains=3, n.adapt=500)
  update(mod, n.iter = iterations)
  out <- coda.samples(mod, c("beta", "intG", "intYr", "interceptMu", "temp1", "temp2", "rain1", "rain2"),
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
  
  saveRDS(outC, file = paste(doSpp, "_GrowthParamsMCMC.rds", sep=""))
  write.csv(gelmDiag, file=paste(doSpp, "_growthGelman.csv", sep=""))
}





