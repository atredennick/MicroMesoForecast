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
#### Read in data by species and make one long data frame -------------
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
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRgrowthCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOgrowthCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMgrowthCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEgrowthCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)

# crowd <- rbind(read.csv("BOGRgrowthCrowding.csv")[,2:3], 
#            read.csv("HECOgrowthCrowding.csv")[,2:3],
#            read.csv("PASMgrowthCrowding.csv")[,2:3],
#            read.csv("POSEgrowthCrowding.csv")[,2:3])
# crowd <- read.csv("HECOgrowthCrowding.csv")[,2:3] #ignore first column of rownames

# Merge crowding and growth data
growD <- merge(growD, crowd, by=c("species", "X"))

#try glm
# fit final mixed effect model: based on what?
library(lme4)
outlm=lmer(log(area.t1)~log(area.t0)+W+pptLag+ppt1+TmeanSpr1+ 
           ppt2+TmeanSpr2+
           ppt1:TmeanSpr1+ppt2:TmeanSpr2+
           (1|Group)+(log(area.t0)|year),data=subset(growD, species=="HECO")) 
summary(outlm)

modelstring="
  model{
  #process model and likelihood
  for(i in 1:nObs){
    mu[i] <- intYr[spp[i],yrs[i]] + intG[spp[i],grp[i]] + beta[spp[i],yrs[i]]*X[i] + 
    nb[spp[i]]*crowd[i] + temp1[spp[i]]*TmeanSpr1[i] + 
    temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + 
    rain2[spp[i]]*ppt2[i] + rainlag[spp[i]]*pptlag[i]+
    rain1Tmn1[spp[i]]*ppt1Tmn1[i] + rain2Tmn2[spp[i]]*ppt2Tmn2[i]
    tau2[i] <- 1/(tau[spp[i]]*exp(tauSize[spp[i]]*mu[i])) 
    tau3[i] <- max(tau2[i],0.00000001)  
    Y[i] ~ dnorm(mu[i], tau3[spp[i]])
  }

  #priors
  for(s in 1:nSpp){
    tau[s] ~ dnorm(0,0.001)
    tauSize[s] ~ dnorm(0,0.001)
    betaSpp[s] ~ dnorm(0,0.001)
    nb[s] ~ dnorm(0, 0.001)
    temp1[s] ~ dnorm(temp1Mu, temp1Var)
    temp2[s] ~ dnorm(temp2Mu, temp2Var)
    rain1[s] ~ dnorm(rain1Mu, rain1Var)
    rain2[s] ~ dnorm(rain2Mu, rain2Var)
    rainlag[s] ~ dnorm(rainlagMu, rainlagVar)
    rain1Tmn1[s] ~ dnorm(rain1Tmn1Mu, rain1Tmn1Var)
    rain2Tmn2[s] ~ dnorm(rain2Tmn2Mu, rain2Tmn2Var)
    intercept[s] ~ dnorm(0, 0.001)
    intVaryY[s] ~ dgamma(0.5, 0.5)
    betaVar[s] ~ dgamma(0.5, 0.5)
    intVarG[s] ~ dgamma(0.001, 0.001) 
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
    }
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }

  temp1Mu ~ dnorm(0,1e-6)
  temp2Mu ~ dnorm(0,1e-6)
  rain1Mu ~ dnorm(0,1e-6)
  rain2Mu ~ dnorm(0,1e-6)
  rainlagMu ~ dnorm(0,1e-6)
  rain1Tmn1Mu ~ dnorm(0,1e-6)
  rain2Tmn2Mu ~ dnorm(0,1e-6)
  temp1Var ~ dgamma(0.5, 0.5)
  temp2Var ~ dgamma(0.5, 0.5)
  rain1Var ~ dgamma(0.5, 0.5)
  rain2Var ~ dgamma(0.5, 0.5)
  rainlagVar ~ dgamma(0.5, 0.5)
  rain1Tmn1Var ~ dgamma(0.5, 0.5)
  rain2Tmn2Var ~ dgamma(0.5, 0.5)
}#end model
"


####
#### Set up data for JAGS model ---------------------------
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
              crowd = growD$W,
              TmeanSpr1 = growD$TmeanSpr1,
              TmeanSpr2 = growD$TmeanSpr2,
              ppt1 = growD$ppt1,
              ppt2 = growD$ppt2,
              pptlag = growD$pptLag,
              ppt1Tmn1 = growD$ppt1*growD$TmeanSpr1,
              ppt2Tmn2 = growD$ppt2*growD$TmeanSpr2)

# Set up some initial values for the MCMC chains (3 chains)
nyrs <- length(unique(growD$year))
nspp <- 4
ngrp <- length(unique(growD$Group))
inits=list(1)
inits[[1]]=list(intercept=rep(1,nspp), intYr=matrix(1, ncol=nyrs, nrow=nspp),
                intG=matrix(0.1,ncol=ngrp, nrow=nspp), betaSpp=rep(1,nspp), betaVar=rep(1,nspp),
                beta=matrix(1,ncol=nyrs,nrow=nspp), nb=rep(-0.2,nspp), intVarG=rep(2,nspp),
                tau=rep(1,nspp), tauSize=rep(1,nspp), intVaryY=rep(1,nspp),
                temp1Mu=1, temp2Mu=1, rain1Mu=1, rain2Mu=1, rainlagMu=1, rain1Tmn1=1, rain2Tmn2=1,
                temp1Var=1, temp2Var=1, rain1Var=1, rain2Var=1, rainlagVar=1,
                temp1=rep(1,nspp), temp2=rep(1,nspp), rain1=rep(1,nspp), rain2=rep(1,nspp), rainlag=rep(1,nspp))
inits[[1]] <- list(tau=rep(1,nspp), tauSize=rep(1,nspp))
inits[[2]]=list(intercept=rep(-1,nspp), intYr=matrix(-1, ncol=nyrs, nrow=nspp),
                intG=matrix(-0.1,ncol=ngrp, nrow=nspp), betaSpp=rep(-1,nspp), betaVar=rep(0.5,nspp),
                beta=matrix(-1,ncol=nyrs,nrow=nspp), nb=rep(0.2,nspp), intVarG=rep(0.5,nspp),
                tau=rep(0.5,nspp), tauSize=rep(0.5,nspp), intVaryY=rep(0.5,nspp),
                temp1Mu=2, temp2Mu=2, rain1Mu=2, rain2Mu=2, rainlagMu=2, rain1Tmn1=2, rain2Tmn2=2,
                temp1Var=2, temp2Var=2, rain1Var=2, rain2Var=2, rainlagVar=1,
                temp1=rep(2,nspp), temp2=rep(2,nspp), rain1=rep(2,nspp), rain2=rep(2,nspp), rainlag=rep(2,nspp))
inits[[2]] <- list(tau=rep(0.5,nspp), tauSize=rep(0.5,nspp))
inits[[3]]=list(intercept=rep(0.1,nspp), intYr=matrix(0.1, ncol=nyrs, nrow=nspp),
                intG=matrix(0.5,ncol=ngrp, nrow=nspp), betaSpp=rep(0.1,nspp), betaVar=rep(0.1,nspp),
                beta=matrix(0.1,ncol=nyrs,nrow=nspp), nb=rep(-0.5,nspp), intVarG=rep(0.2,nspp),
                tau=rep(0.1,nspp), tauSize=rep(0.1,nspp), intVaryY=rep(0.1,nspp),
                temp1Mu=0.1, temp2Mu=0.1, rain1Mu=0.1, rain2Mu=0.1, rainlagMu=0.1, rain1Tmn1=0.1, rain2Tmn2=0.1,
                temp1Var=0.1, temp2Var=0.1, rain1Var=0.1, rain2Var=0.1, rainlagVar=0.1,
                temp1=rep(0.1,nspp), temp2=rep(0.1,nspp), rain1=rep(0.1,nspp), rain2=rep(0.1,nspp), rainlag=rep(0.1,nspp))
inits[[3]] <- list(tau=rep(0.1,nspp), tauSize=rep(0.1,nspp))

####
#### Run MCMC from JAGS ------------------------
####
iterations <- 2500
adapt <- 500

mod <- jags.model(textConnection(modelstring), data=dataJ, n.chains=length(inits), 
                  n.adapt=adapt, inits=inits)
update(mod, n.iter = (iterations))
# out <- coda.samples(mod, c("intYr", "beta", "intG", "nb", "temp1", "temp2", 
#                            "rain1", "rain2", "intercept", "betaSpp",
#                            "rain1Tmn1", "rain2Tmn2", "tau", "tauSize"),
#                     n.iter=iterations, n.thin=1)
out <- coda.samples(mod, c( "nb", "temp1"),
                    n.iter=iterations, n.thin=5)
plot(out)
# dic <- jags.samples(mod, c("deviance"),
#                     n.iter=iterations, n.thin=10)

####
#### Check for convergence ---------------------
####
gelmDiag <- gelman.diag(out)
# heidel.diag(out)
# gelman.plot(out)

# pdf("growthOutPlots.pdf")
# plot(out, auto.layout=FALSE)
# dev.off()

####
#### Convert to dataframe for export and get other summaries ------------------
####
outC <- rbind(out[[1]], 
              out[[2]], 
              out[[3]])

# head(outC)
outD <- as.data.frame(out[[1]])
plot(outD$temp1, type="l")

# outC <- rbind(out[[1]][(iterations-999):iterations,], 
#               out[[2]][(iterations-999):iterations,], 
#               out[[3]][(iterations-999):iterations,])

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
