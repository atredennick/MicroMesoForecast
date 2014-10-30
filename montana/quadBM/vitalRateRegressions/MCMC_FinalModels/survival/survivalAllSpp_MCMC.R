#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)

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
  
  #Survival observations
  survD <- subset(sppD,lag.cover>0)
  survD$survives <- ifelse(survD$totCover>0,1,0)
  survD$yearID <- survD$year #for random year offset on intercept
  survD$group <- substring(survD$quad, 1, 1)
  survD$percCover <- survD$totCover/10000
  survD$percLagCover <- survD$lag.cover/10000
  
  
  ####
  #### Set up data structure for JAGS
  ####
  nGrp <- length(unique(survD$group))
  nYrs <- length(unique(survD$year))
  nObs <- nrow(survD)
  size <- survD$percLagCover
  S <- survD$survives
  yrs <- (survD$year)-32
  grp <- as.numeric(as.factor(survD$group))
  TmeanSpr1 <- survD$TmeanSpr1
  TmeanSpr2 <- survD$TmeanSpr2
  ppt1 <- survD$ppt1
  ppt2 <- survD$ppt2
  
  ####
  #### Run MCMC
  ####
  iterations <- 30000
  adapt <- 10000
  dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, size=size, S=S, yrs=yrs, grp=grp,
                TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2)
  mod <- jags.model("survivalAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
  update(mod, n.iter = (iterations/2))
  out <- coda.samples(mod, c("beta", "intG", "intYr", "interceptMu", "temp1", "temp2", "rain1", "rain2"),
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





